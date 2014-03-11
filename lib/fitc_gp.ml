(* File: fitc_gp.ml

   Sized GPR - OCaml-GPR with static size checking of operations on matrices

   [Authors of Sized GPR]
     Copyright (C) 2014-  Akinori ABE
     email: abe@kb.ecei.tohoku.ac.jp

   [Authors of OCaml-GPR]
     Copyright (C) 2009-  Markus Mottl
     email: markus.mottl@gmail.com
     WWW:   http://www.ocaml.info

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*)

open Gpr_utils
open Interfaces

open Core.Std
open Bigarray
open Slap.D (*! RID *)

module type Sig = functor (Spec : Specs.Eval) ->
  Sigs.Eval with module Spec = Spec

(* Computations shared by FIC and FITC, and standard and variational
   version *)
module Make_common (Spec : Specs.Eval) = struct
  module Spec = Spec

  open Spec

  let jitter = !cholesky_jitter

  (* Evaluation of inducing points *)
  module Inducing = struct
    type ('D, 'd, 'm) t = { (*! ITP *)
      kernel : ('D, 'd, 'm) Kernel.t; (*! ITP *)
      points : ('d, 'm) Spec.Inducing.t; (*! ITP *)
      km : ('m, 'm, Slap.cnt) mat; (*! ITP *)
      chol_km : ('m, 'm, Slap.cnt) mat; (*! ITP *)
      log_det_km : float;
    }

    let check_n_inducing ~n_inducing inputs =
      let n_inputs = Slap.Size.to_int (Spec.Inputs.get_n_points inputs) (*! S2I *) in
      if n_inputs < 1 || (Slap.Size.to_int n_inducing) (*! S2I *) > n_inputs then
        failwithf
          "Gpr.Fitc_gp.Make_common.check_n_inducing: \
          violating 1 <= n_inducing (%d) <= n_inputs (%d)"
          (Slap.Size.to_int n_inducing) (*! S2I *) n_inputs ()

    let calc_internal kernel points km =
      let chol_km = lacpy ~uplo:`U km in
      potrf ~jitter chol_km;
      { kernel; points; km; chol_km; log_det_km = log_det chol_km }

    let calc kernel points =
      calc_internal kernel points (Spec.Inducing.calc_upper kernel points)

    let choose kernel inputs indexes =
      let chosen_inputs = Spec.Inputs.choose_subset inputs indexes in
      Spec.Inputs.create_inducing kernel chosen_inputs

    let choose_n_first_inputs kernel inputs ~n_inducing =
      check_n_inducing ~n_inducing inputs;
      let indexes = Int_vec.create n_inducing in
      for i = 1 to Slap.Size.to_int n_inducing (*! S2I *) do Slap.Vec.set_dyn indexes i i (*! IDX *) done;
      choose kernel inputs indexes

    let choose_n_random_inputs
          ?(rnd_state = Random.State.default) kernel inputs ~n_inducing =
      check_n_inducing ~n_inducing inputs;
      let n_inputs = Spec.Inputs.get_n_points inputs in
      let indexes = Int_vec.create n_inputs in
      for i = 1 to Slap.Size.to_int n_inputs (*! S2I *) do Slap.Vec.set_dyn indexes i i (*! IDX *)  done;
      for i = 1 to Slap.Size.to_int n_inducing (*! S2I *)  do
        let rnd_index = Random.State.int rnd_state ((Slap.Size.to_int n_inputs) (*! S2I *) - i + 1) + 1 in
        let tmp = Slap.Vec.get_dyn indexes rnd_index (*! IDX *) in
        Slap.Vec.set_dyn indexes rnd_index (Slap.Vec.get_dyn indexes i); (*! IDX *)
        Slap.Vec.set_dyn indexes i tmp; (*! IDX *)
      done;
      let indexes = Slap.Vec.subvec_dyn ~ofsx:1 n_inducing indexes (*! RID *) in
      choose kernel inputs indexes

    let get_kernel inducing = inducing.kernel
    let get_points inducing = inducing.points
  end

  (* Evaluation of one input point *)
  module Input = struct
    type ('D, 'd, 'm) t = { inducing : ('D, 'd, 'm) Inducing.t; point : 'D Spec.Input.t; k_m : ('m, Slap.cnt) vec } (*! ITP *)

    let calc inducing point =
      let { Inducing.kernel; points = inducing_points } = inducing in
      { inducing; point; k_m = Spec.Input.eval kernel point inducing_points }
  end

  (* Evaluation of input points *)
  module Inputs = struct
    type ('D, 'd, 'm, 'n) t = { inducing : ('D, 'd, 'm) Inducing.t; points : ('D, 'n) Inputs.t; knm : ('n, 'm, Slap.cnt) mat } (*! ITP *)

    let calc_internal points inducing knm = { inducing; points; knm }

    let calc points inducing =
      let { Inducing.kernel; points = inducing_points } = inducing in
      let knm =
        Inputs.calc_cross kernel ~inputs:points ~inducing:inducing_points
      in
      { inducing; points; knm }

    let get_kernel t = t.inducing.Inducing.kernel

    let calc_diag inputs = Inputs.calc_diag (get_kernel inputs) inputs.points
    let calc_upper inputs = Inputs.calc_upper (get_kernel inputs) inputs.points

    let create_default_kernel inputs ~n_inducing =
      Kernel.create (
        Spec.Inputs.create_default_kernel_params inputs ~n_inducing)

    let get_chol_km t = t.inducing.Inducing.chol_km
    let get_log_det_km t = t.inducing.Inducing.log_det_km
    let get_knm t = t.knm
    let get_points t = t.points
  end

  (* Model computations shared by standard and variational version *)
  module Common_model = struct
    type ('D, 'd, 'm, 'n) t = { (*! ITP *)
      sigma2 : float;
      inputs : ('D, 'd, 'm, 'n) Inputs.t; (*! ITP *)
      kn_diag : ('n, Slap.cnt) vec; (*! ITP *)
      v_mat : ('n, 'm, Slap.cnt) mat; (*! ITP *)
      r_vec : ('n, Slap.cnt) vec; (*! ITP *)
      is_vec : ('n, Slap.cnt) vec; (*! ITP *)
      sqrt_is_vec : ('n, Slap.cnt) vec; (*! ITP *)
      q_mat : (('n, 'm) Slap.Size.add, 'm, Slap.cnt) mat; (*! ITP *)
      r_mat : ('m, 'm, Slap.cnt) mat; (*! ITP *)
      l1 : float;
    }

    type 'n co_variance_coeffs = ('n, 'n, Slap.cnt) mat * ('n, 'n, Slap.cnt) mat (*! ITP *)

    let check_sigma2 sigma2 =
      if sigma2 < 0. then failwith "Model.check_sigma2: sigma2 < 0"

    let calc_internal inputs sigma2 ~kn_diag ~v_mat ~r_vec =
      check_sigma2 sigma2;
      let n = Mat.dim1 v_mat in
      let m = Mat.dim2 v_mat in
      let is_vec = Vec.create n in
      let log_det_s_vec =
        let rec loop log_det_s_vec i =
          if i = 0 then log_det_s_vec
          else
            let s_vec_i = (Vec.get_dyn r_vec i) (*! IDX *) +. sigma2 in
            let is_vec_i = 1. /. s_vec_i in
            Vec.set_dyn is_vec i is_vec_i; (*! IDX *)
            loop (log_det_s_vec +. log s_vec_i) (i - 1)
        in
        loop 0. (Slap.Size.to_int n) (*! S2I *)
      in
      let sqrt_is_vec = Vec.sqrt is_vec in
      let nm = (Slap.Size.to_int n) + (Slap.Size.to_int m) (*! S2I *) in
      let n1 = (Slap.Size.to_int n) (*! S2I *) + 1 in
      let q_mat = Mat.create (Slap.Size.add n m) (*! SOP *) m in
      for c = 1 to Slap.Size.to_int m (*! S2I *) do
        for r = n1 + c to nm do Mat.set_dyn q_mat r c 0. (*! IDX *) done;
      done;
      let t = Inputs.get_knm inputs in (*! SUB *)
      ignore (lacpy t ~b:(Slap.Mat.submat_dyn (Mat.dim1 t) (Mat.dim2 t) q_mat)); (*! SUB *)
      Mat.scal_rows sqrt_is_vec (Slap.Mat.submat_dyn n m q_mat); (*! SUB *)
      let chol_km = Inputs.get_chol_km inputs in
      ignore (lacpy ~uplo:`U chol_km ~b:(Slap.Mat.submat_dyn (Mat.dim1 chol_km) (Mat.dim2 chol_km) ~ar:n1 q_mat)); (*! SUB *)
      let tau = geqrf q_mat in
      let r_mat = lacpy ~uplo:`U (Slap.Mat.submat_dyn m m q_mat) in (*! SUB *)
      orgqr_dyn ~tau q_mat; (*! RID *)
      let log_det_r =
        let rec loop log_det_r r =
          if r = 0 then log_det_r +. log_det_r
          else
            let el =
              let el = Mat.get_dyn r_mat r r (*! IDX *) in
              if el > 0. then el
              else
                (* Cannot happen with LAPACK version 3.2 and greater *)
                let neg_el = -. el in
                Mat.set_dyn r_mat r r neg_el; (*! IDX *)
                for c = r + 1 to Slap.Size.to_int m (*! S2I *) do Mat.set_dyn r_mat r c (-. (Mat.get_dyn r_mat r c)) (*! IDX *) done;
                Mat.scal ~-.1. (Slap.Mat.submat_dyn n Slap.Size.one (*! SC *) ~ac:r q_mat); (*! SUB *)
                neg_el
            in
            loop (log_det_r +. log el) (r - 1)
        in
        loop 0. (Slap.Size.to_int m) (*! S2I *)
      in
      let l1 =
        let log_det_km = Inputs.get_log_det_km inputs in
        -0.5 *. (log_det_r -. log_det_km +. log_det_s_vec +. float (Slap.Size.to_int n) (*! S2I *) *. log_2pi)
      in
      {
        inputs; sigma2; kn_diag; v_mat; r_vec;
        is_vec; sqrt_is_vec; q_mat; r_mat; l1;
      }

    let calc_r_vec ~kn_diag ~v_mat =
      Mat.syrk_diag ~trans:Slap.Common.normal (*! IF *) ~alpha:~-.1. v_mat ~beta:1. ~y:(copy kn_diag)

    let calc_with_kn_diag inputs sigma2 kn_diag =
      let v_mat = lacpy inputs.Inputs.knm in
      trsm ~side:Slap.Common.right (*! RF *) ~transa:Slap.Common.normal (*! IF *) v_mat ~a:(Inputs.get_chol_km inputs);
      let r_vec = calc_r_vec ~kn_diag ~v_mat in
      calc_internal inputs sigma2 ~kn_diag ~v_mat ~r_vec

    let calc inputs ~sigma2 =
      calc_with_kn_diag inputs sigma2 (Inputs.calc_diag inputs)

    let update_sigma2 { inputs; kn_diag; v_mat; r_vec } sigma2 =
      check_sigma2 sigma2;
      calc_internal inputs sigma2 ~kn_diag ~v_mat ~r_vec

    let calc_log_evidence model = model.l1

    let get_v_mat model = model.v_mat
    let get_inducing model = model.inputs.Inputs.inducing
    let get_inducing_points model = Inducing.get_points (get_inducing model)
    let get_inputs model = model.inputs
    let get_input_points model = model.inputs.Inputs.points
    let get_kernel model = (get_inducing model).Inducing.kernel
    let get_sigma2 model = model.sigma2
    let get_chol_km model = Inputs.get_chol_km model.inputs
    let get_r_vec model = model.r_vec
    let get_is_vec model = model.is_vec
    let get_sqrt_is_vec model = model.sqrt_is_vec
    let get_km model = (get_inducing model).Inducing.km
    let get_knm model = model.inputs.Inputs.knm
    let get_kn_diag model = model.kn_diag
    let get_q_mat model = model.q_mat
    let get_r_mat model = model.r_mat

    let calc_co_variance_coeffs model = get_chol_km model, model.r_mat
  end

  (* Model computation (variational version) *)
  module Variational_model = struct
    include Common_model

    let from_common ({ r_vec; is_vec; l1 } as model) =
      { model with l1 = l1 +. -0.5 *. dot ~x:is_vec r_vec }

    let calc_with_kn_diag inputs sigma2 kn_diag =
      from_common (calc_with_kn_diag inputs sigma2 kn_diag)

    let calc inputs ~sigma2 = from_common (calc inputs ~sigma2)
    let update_sigma2 model sigma2 = from_common (update_sigma2 model sigma2)
  end

  (* Trained models *)
  module Trained = struct
    type ('D, 'd, 'm, 'n) t = { (*! ITP *)
      model : ('D, 'd, 'm, 'n) Common_model.t; (*! ITP *)
      y : ('n, Slap.cnt) vec; (*! ITP *)
      coeffs : ('m, Slap.cnt) vec; (*! ITP *)
      l : float;
    }

    let calc_internal model ~y ~coeffs ~l2 =
      { model; y; coeffs; l = model.Common_model.l1 +. l2 }

    let prepare_internal model ~y =
      let sqrt_is_vec = model.Common_model.sqrt_is_vec in
      let n = Vec.dim sqrt_is_vec in
      let y_ = Vec.mul y sqrt_is_vec in
      y_, gemv ~trans:Slap.Common.trans (*! RF *) (Slap.Mat.submat_dyn n (Mat.dim2 model.Common_model.q_mat) model.Common_model.q_mat) y_ (*! SUB *)

    let calc model ~targets:y =
      let y_, qt_y_ = prepare_internal model ~y in
      let l2 = -0.5 *. (Vec.sqr_nrm2 y_ -. Vec.sqr_nrm2 qt_y_) in
      trsv ~trans:Slap.Common.normal (*! IF *) model.Common_model.r_mat qt_y_;
      calc_internal model ~y ~coeffs:qt_y_ ~l2

    let calc_mean_coeffs trained = trained.coeffs
    let calc_log_evidence trained = trained.l

    let calc_means trained =
      gemv ~trans:Slap.Common.normal (*! IF *) (Common_model.get_knm trained.model) trained.coeffs

    let get_inducing trained = Common_model.get_inducing trained.model
    let get_targets (trained : (_, _, _, _) t) (*! ITP *) = trained.y
    let get_model trained = trained.model
  end

  module Stats = struct
      type t = {
        n_samples : int;
        target_variance : float;
        sse : float;
        mse : float;
        rmse : float;
        smse : float;
        msll : float;
        mad : float;
        maxad : float;
      }

    let calc_n_samples { Trained.y } = Vec.dim y

    let calc_target_variance { Trained.y } =
      Vec.sqr_nrm2 y /. float (Slap.Size.to_int (Vec.dim y)) (*! S2I *)

    let calc_sse trained =
      let means = Trained.calc_means trained in
      Vec.ssqr_diff trained.Trained.y means

    let calc_mse trained = calc_sse trained /. float (Slap.Size.to_int (calc_n_samples trained)) (*! S2I *)
    let calc_rmse trained = sqrt (calc_mse trained)
    let calc_smse trained = calc_mse trained /. calc_target_variance trained

    let calc_prior_l target_variance =
      -0.5 *. log (2. *. pi *. target_variance) -. 0.5

    let calc_msll trained =
      let prior_l = calc_prior_l (calc_target_variance trained) in
      prior_l -. trained.Trained.l /. float (Slap.Size.to_int (calc_n_samples trained)) (*! S2I *)

    let calc_mad ({ Trained.y } as trained) =
      let n_samples = Slap.Size.to_int (Vec.dim y) (*! S2I *) in
      let f_samples = float n_samples in
      let means = Trained.calc_means trained in
      let rec loop madsum i =
        if i = 0 then madsum /. f_samples
        else loop (madsum +. Float.abs ((Vec.get_dyn y i) -. (Vec.get_dyn means i))) (*! IDX *) (i - 1)
      in
      loop 0. n_samples

    let calc_maxad ({ Trained.y } as trained) =
      let means = Trained.calc_means trained in
      let rec loop maxad i =
        if i = 0 then maxad
        else loop (max maxad (Float.abs ((Vec.get_dyn y i) -. (Vec.get_dyn means i)))) (*! IDX *) (i - 1)
      in
      loop 0. (Slap.Size.to_int (Vec.dim y)) (*! S2I *)

    let calc ({ Trained.y; l } as trained) =
      let n_samples = Slap.Size.to_int (Vec.dim y) (*! S2I *) in
      let f_samples = float n_samples in
      let target_variance = calc_target_variance trained in
      let means = Trained.calc_means trained in
      let sse = Vec.ssqr_diff y means in
      let mse = sse /. f_samples in
      let rmse = sqrt mse in
      let smse = mse /. target_variance in
      let prior_l = calc_prior_l target_variance in
      let msll = prior_l -. l /. f_samples in
      let mad, maxad =
        let rec loop ~madsum ~maxad i =
          if i = 0 then madsum /. f_samples, maxad
          else
            let ad = Float.abs ((Vec.get_dyn y i) -. (Vec.get_dyn means i)) (*! IDX *) in
            loop ~madsum:(madsum +. ad) ~maxad:(max maxad ad) (i - 1)
        in
        loop ~madsum:0. ~maxad:0. n_samples
      in
      { n_samples; target_variance; sse; mse; rmse; smse; msll; mad; maxad }
  end

  module Mean_predictor = struct
    type ('m, 'n) t = { inducing : ('m, 'n) Spec.Inducing.t; coeffs : ('n, Slap.cnt) vec } (*! ITP *)

    let calc_trained trained =
      {
        inducing = Inducing.get_points (Trained.get_inducing trained);
        coeffs = Trained.calc_mean_coeffs trained;
      }

    let calc inducing ~coeffs =
      (* if Spec.Inducing.get_n_points inducing <> Vec.dim coeffs then *) (*! RMDC *)
      (*   failwith *) (*! RMDC *)
      (*     "Mean_predictor.calc: number of inducing points disagrees with \ *) (*! RMDC *)
      (*     dimension of coefficients" *) (*! RMDC *)
      (* else *) { inducing; coeffs } (*! RMDC *)

    let get_inducing t = t.inducing
    let get_coeffs t = t.coeffs
  end

  (* Prediction of mean for one input point *)
  module Mean = struct
    type 'n t = { point : 'n Spec.Input.t; value : float } (*! ITP *)

    let calc mean_predictor { Input.inducing = input_inducing; k_m; point } =
      if
        not (
          phys_equal
            mean_predictor.Mean_predictor.inducing
            (Inducing.get_points input_inducing))
      then
        failwith
          "Mean.calc: mean predictor and input disagree about inducing points"
      else { point; value = dot ~x:k_m mean_predictor.Mean_predictor.coeffs }

    let get mean = mean.value
  end

  (* Prediction of means for several input points *)
  module Means = struct
    type ('m, 'n) t = { points : ('m, 'n) Spec.Inputs.t; values : ('n, Slap.cnt) vec } (*! ITP *)

    let calc mean_predictor { Inputs.points; inducing; knm } =
      if
        not (
          phys_equal
            mean_predictor.Mean_predictor.inducing
            (Inducing.get_points inducing))
      then
        failwith "Means.calc: trained and inputs disagree about inducing points"
      else { points; values = gemv ~trans:Slap.Common.normal (*! IF *) knm mean_predictor.Mean_predictor.coeffs }

    let get means = means.values
  end

  module Co_variance_predictor = struct
    type ('D, 'd, 'm, 'n) t = { (*! ITP *)
      kernel : ('D, 'd, 'm) Spec.Kernel.t; (*! ITP *)
      inducing : ('n, 'm) Spec.Inducing.t; (*! ITP *)
      chol_km : ('m, 'm, Slap.cnt) mat; (*! ITP *)
      r_mat : ('m, 'm, Slap.cnt) mat; (*! ITP *)
    }

    let calc_model model =
      {
        kernel = Common_model.get_kernel model;
        inducing = Inducing.get_points (Common_model.get_inducing model);
        chol_km = Common_model.get_chol_km model;
        r_mat = Common_model.get_r_mat model;
      }

    let calc kernel inducing (chol_km, r_mat) =
      { kernel; inducing; chol_km; r_mat }
  end

  (* Prediction of variance for one input point *)
  module Variance = struct
    type 'n t = { point : 'n Spec.Input.t; variance : float; sigma2 : float } (*! ITP *)

    let calc co_variance_predictor ~sigma2 { Input.inducing; point; k_m } =
      if
        not (
          phys_equal
            co_variance_predictor.Co_variance_predictor.inducing
            (Inducing.get_points inducing))
      then
        failwith
          "Variance.calc: \
          co-variance predictor and input disagree about inducing points"
      else
        let variance =
          let
            { Co_variance_predictor.kernel; chol_km; r_mat } =
              co_variance_predictor
          in
          let tmp = copy k_m in
          trsv ~trans:Slap.Common.trans (*! RF *) chol_km tmp;
          let k = Vec.sqr_nrm2 tmp in
          let tmp = copy k_m ~y:tmp in
          trsv ~trans:Slap.Common.trans (*! RF *) r_mat tmp;
          let b = Vec.sqr_nrm2 tmp in
          let prior_variance = Spec.Input.eval_one kernel point in
          prior_variance -. (k -. b)
        in
        { point; variance; sigma2 }

    let get ?predictive t =
      match predictive with
      | None | Some true -> t.variance +. t.sigma2
      | Some false -> t.variance
  end

  (* Prediction of variance for several input points *)
  module Variances = struct
    type ('m, 'n) t = { points : ('m, 'n) Spec.Inputs.t; variances : ('n, Slap.cnt) vec; sigma2 : float } (*! ITP *)

    let calc_model_inputs model =
      let variances =
        let tmp = lacpy (Common_model.get_knm model) in
        trsm ~side:Slap.Common.right (*! RF *) ~transa:Slap.Common.normal (*! IF *) tmp ~a:model.Common_model.r_mat;
        Mat.syrk_diag ~trans:Slap.Common.normal (*! IF *) tmp ~beta:1. ~y:(copy (Common_model.get_r_vec model))
      in
      let sigma2 = Common_model.get_sigma2 model in
      { points = Common_model.get_input_points model; variances; sigma2 }

    let calc cvp ~sigma2 inputs =
      if
        not (
          phys_equal
            cvp.Co_variance_predictor.inducing
            (Inducing.get_points inputs.Inputs.inducing))
      then
        failwith
          "Variances.calc: \
          co-variance predictor and inputs disagree about inducing points"
      else
        let { Inputs.points; knm = ktm } = inputs in
        let variances =
          let y = Inputs.calc_diag inputs in
          let tmp = lacpy ktm in
          trsm ~side:Slap.Common.right (*! RF *) ~transa:Slap.Common.normal (*! IF *) tmp ~a:cvp.Co_variance_predictor.chol_km;
          let y = Mat.syrk_diag ~trans:Slap.Common.normal (*! IF *) ~alpha:~-.1. tmp ~beta:1. ~y in
          let tmp = lacpy ktm ~b:tmp in
          trsm ~side:Slap.Common.right (*! RF *) ~transa:Slap.Common.normal (*! IF *) tmp ~a:cvp.Co_variance_predictor.r_mat;
          Mat.syrk_diag ~trans:Slap.Common.normal (*! IF *) tmp ~beta:1. ~y
        in
        { points; variances; sigma2 }

    let get_common ?predictive ~variances ~sigma2 =
      match predictive with
      | None | Some true ->
          let predictive_variances = Vec.make (Vec.dim variances) sigma2 in
          axpy ~x:variances predictive_variances;
          predictive_variances
      | Some false -> variances

    let get ?predictive { variances; sigma2 } =
      get_common ?predictive ~variances ~sigma2
  end

  (* Computations for predicting covariances shared by FIC and
     FITC, and standard and variational version *)
  module Common_covariances = struct
    type ('m, 'n) t = { points : ('m, 'n) Spec.Inputs.t; covariances : ('n, 'n, Slap.cnt) mat; sigma2 : float } (*! ITP *)

    let check_inducing ~loc co_variance_predictor inputs =
      if
        not (
          phys_equal
            co_variance_predictor.Co_variance_predictor.inducing
            (Inducing.get_points inputs.Inputs.inducing))
      then
        failwithf
          "%s_covariances.calc: \
          co-variance predictor and inputs disagree about inducing points"
          loc ()

    let get_common ?predictive ~covariances ~sigma2 =
      match predictive with
      | None | Some true ->
          let res = lacpy ~uplo:`U covariances in
          for i = 1 to Slap.Size.to_int (Mat.dim1 res) (*! S2I *) do Mat.set_dyn res i i ((Mat.get_dyn res i i) +. sigma2) (*! IDX *) done;
          res
      | Some false -> covariances

    let get ?predictive { covariances; sigma2 } =
      get_common ?predictive ~covariances ~sigma2

    let get_variances { points; covariances; sigma2 } =
      { Variances.points; variances = Mat.copy_diag covariances; sigma2 }
  end

  (* Predicting covariances with FITC (standard or variational) *)
  module FITC_covariances = struct
    include Common_covariances

    let calc_model_inputs model =
      let covariances = Inputs.calc_upper model.Common_model.inputs in
      let v_mat = model.Common_model.v_mat in
      ignore (syrk ~trans:Slap.Common.normal (*! IF *) ~alpha:~-.1. v_mat ~beta:1. ~c:covariances);
      let q_mat = model.Common_model.q_mat in
      let n = Mat.dim1 v_mat in
      ignore (syrk ~trans:Slap.Common.normal (*! IF *) (Slap.Mat.submat_dyn n (Mat.dim2 q_mat) q_mat) ~beta:1. ~c:covariances); (*! SUB *)
      let points = Common_model.get_input_points model in
      let sigma2 = Common_model.get_sigma2 model in
      { points; covariances; sigma2 }

    let calc co_variance_predictor ~sigma2 inputs =
      check_inducing ~loc:"FITC" co_variance_predictor inputs;
      let { Co_variance_predictor.chol_km; r_mat } = co_variance_predictor in
      let covariances = Inputs.calc_upper inputs in
      let { Inputs.points; knm = ktm } = inputs in
      let covariances =
        let tmp = lacpy ktm in
        trsm ~side:Slap.Common.right (*! RF *) ~transa:Slap.Common.normal (*! IF *) tmp ~a:chol_km;
        ignore (syrk ~trans:Slap.Common.normal (*! IF *) ~alpha:~-.1. tmp ~c:covariances);
        let tmp = lacpy ktm ~b:tmp in
        trsm ~side:Slap.Common.right (*! RF *) ~transa:Slap.Common.normal (*! IF *) tmp ~a:r_mat;
        syrk ~trans:Slap.Common.normal (*! IF *) tmp ~c:covariances;
      in
      { points; covariances; sigma2 }
  end

  (* Predicting covariances with FIC (standard or variational) *)
  module FIC_covariances = struct
    include Common_covariances

    let calc_common ~points ~sigma2 ~q_mat ~r_vec =
      let n = Vec.dim r_vec in
      let covariances = syrk ~trans:Slap.Common.normal (*! IF *) (Slap.Mat.submat_dyn n (Mat.dim2 q_mat) q_mat) in (*! SUB *)
      for i = 1 to Slap.Size.to_int n (*! S2I *) do
        Mat.set_dyn covariances i i ((Mat.get_dyn covariances i i) +. (Vec.get_dyn r_vec i)) (*! IDX *)
      done;
      { points; covariances; sigma2 }

    let calc_model_inputs model =
      let r_vec = model.Common_model.r_vec in
      let q_mat = model.Common_model.q_mat in
      let points = Common_model.get_input_points model in
      let sigma2 = Common_model.get_sigma2 model in
      calc_common ~points ~sigma2 ~q_mat ~r_vec

    let calc co_variance_predictor ~sigma2 ({ Inputs.knm = ktm } as inputs) =
      check_inducing ~loc:"FIC" co_variance_predictor inputs;
      let kt_diag = Inputs.calc_diag inputs in
      let r_vec = Mat.syrk_diag ~trans:Slap.Common.normal (*! IF *) ~alpha:~-.1. ktm ~beta:1. ~y:kt_diag in
      let q_mat = lacpy ktm in
      let r_mat = co_variance_predictor.Co_variance_predictor.r_mat in
      trsm ~side:Slap.Common.right (*! RF *) ~transa:Slap.Common.normal (*! IF *) q_mat ~a:r_mat;
      let points = Inputs.get_points inputs in
      calc_common ~points ~sigma2 ~q_mat ~r_vec
  end

  (* Computations for sampling the marginal posterior GP distribution
     shared by standard and variational version *)
  module Common_sampler = struct
    type t = { mean : float; stddev : float }

    let calc ~loc ?predictive mean variance =
      if not (phys_equal mean.Mean.point variance.Variance.point) then
        failwith (
          loc ^ ".Sampler: mean and variance disagree about input point");
      let used_variance =
        match predictive with
        | None | Some true ->
            variance.Variance.variance +. variance.Variance.sigma2
        | Some false -> variance.Variance.variance
      in
      { mean = mean.Mean.value; stddev = sqrt used_variance }

    let sample ?(rng = default_rng) sampler =
      let noise = Gsl.Randist.gaussian_ziggurat rng ~sigma:sampler.stddev in
      sampler.mean +. noise

    let samples ?(rng = default_rng) sampler ~n =
      Vec.init n (fun _ -> sample ~rng sampler)
  end

  (* Computations for sampling the posterior GP distribution shared
     by FIC and FITC, and standard and variational version *)
  module Common_cov_sampler = struct
    type 'n t = { means : ('n, Slap.cnt) vec; cov_chol : ('n, 'n, Slap.cnt) mat } (*! ITP *)

    let calc ~loc ?predictive means covariances =
      let module Covariances = Common_covariances in
      if not (phys_equal means.Means.points covariances.Covariances.points) then
        failwith (
          loc ^
          ".Cov_sampler: means and covariances disagree about input points");
      let cov_chol = lacpy ~uplo:`U covariances.Covariances.covariances in
      begin
        match predictive with
        | None | Some true ->
            let sigma2 = covariances.Covariances.sigma2 in
            for i = 1 to Slap.Size.to_int (Mat.dim1 cov_chol) (*! S2I *) do
              Mat.set_dyn cov_chol i i ((Mat.get_dyn cov_chol i i) +. sigma2) (*! IDX *)
            done
        | Some false -> ()
      end;
      potrf ~jitter cov_chol;
      { means = means.Means.values; cov_chol }

    let sample ?(rng = default_rng) samplers =
      let n = Vec.dim samplers.means in
      let sample =
        Vec.init n (fun _ -> Gsl.Randist.gaussian_ziggurat rng ~sigma:1.)
      in
      trmv ~trans:Slap.Common.trans (*! RF *) samplers.cov_chol sample;
      axpy ~x:samplers.means sample;
      sample

    let samples ?(rng = default_rng) { means; cov_chol } ~n =
      let n_means = Vec.dim means in
      let samples =
        Mat.init_cols n_means n (fun _ _ ->
          Gsl.Randist.gaussian_ziggurat rng ~sigma:1.)
      in
      trmm ~side:Slap.Common.left (*! IF *) ~transa:Slap.Common.trans (*! RF *) ~a:cov_chol samples;
      for col = 1 to Slap.Size.to_int n (*! S2I *) do
        for row = 1 to Slap.Size.to_int n_means (*! S2I *) do
          Mat.set_dyn samples row col ((Mat.get_dyn samples row col) +. (Vec.get_dyn means row)) (*! IDX *)
        done
      done;
      samples
  end
end

let fitc_loc = "FITC"
let fic_loc = "FIC"
let variational_fitc_loc = "Variational_FITC"
let variational_fic_loc = "Variational_FIC"

module Make_FITC (Spec : Specs.Eval) = struct
  include Make_common (Spec)

  module Model = Common_model
  module Covariances = FITC_covariances

  module Sampler = struct
    include Common_sampler
    let calc ?predictive mean variance = calc ~loc:fitc_loc ?predictive mean variance (*! ETA *)
  end

  module Cov_sampler = struct
    include Common_cov_sampler
    let calc ?predictive means covariances = calc ~loc:fitc_loc ?predictive means covariances (*! ETA *)
  end
end

module Make_FIC (Spec : Specs.Eval) = struct
  include Make_common (Spec)

  module Model = Common_model
  module Covariances = FIC_covariances

  module Sampler = struct
    include Common_sampler
    let calc ?predictive mean variance = calc ~loc:fic_loc ?predictive mean variance (*! ETA *)
  end

  module Cov_sampler = struct
    include Common_cov_sampler
    let calc ?predictive means covariances = calc ~loc:fic_loc ?predictive means covariances (*! ETA *)
  end
end

module Make_variational_FITC (Spec : Specs.Eval) = struct
  include Make_common (Spec)

  module Model = Variational_model
  module Covariances = FITC_covariances

  module Sampler = struct
    include Common_sampler
    let calc ?predictive mean variance = calc ~loc:variational_fitc_loc ?predictive mean variance (*! ETA *)
  end

  module Cov_sampler = struct
    include Common_cov_sampler
    let calc ?predictive means covariances = calc ~loc:variational_fitc_loc ?predictive means covariances (*! ETA *)
  end
end

module Make_variational_FIC (Spec : Specs.Eval) = struct
  include Make_common (Spec)

  module Model = Variational_model
  module Covariances = FIC_covariances

  module Sampler = struct
    include Common_sampler
    let calc ?predictive mean variance = calc ~loc:variational_fic_loc ?predictive mean variance (*! ETA *)
  end

  module Cov_sampler = struct
    include Common_cov_sampler
    let calc ?predictive means covariances = calc ~loc:variational_fic_loc ?predictive means covariances (*! ETA *)
  end
end

module Make (Spec : Specs.Eval) = struct
  module type Sig = Sigs.Eval with module Spec = Spec

  module Common = Make_common (Spec)

  module FITC = struct
    include Common
    module Model = Common_model
    module Covariances = FITC_covariances

    module Sampler = struct
      include Common_sampler
      let calc ?predictive mean variance = calc ~loc:fitc_loc ?predictive mean variance (*! ETA *)
    end

    module Cov_sampler = struct
      include Common_cov_sampler
      let calc ?predictive means covariances = calc ~loc:fitc_loc ?predictive means covariances (*! ETA *)
    end
  end

  module FIC = struct
    include Common
    module Model = Common_model
    module Covariances = FIC_covariances

    module Sampler = struct
      include Common_sampler
      let calc ?predictive mean variance = calc ~loc:fic_loc ?predictive mean variance (*! ETA *)
    end

    module Cov_sampler = struct
      include Common_cov_sampler
      let calc ?predictive means covariances = calc ~loc:fic_loc ?predictive means covariances (*! ETA *)
    end
  end

  module Variational_FITC = struct
    include Common
    module Model = Variational_model
    module Covariances = FITC_covariances

    module Sampler = struct
      include Common_sampler
      let calc ?predictive mean variance = calc ~loc:variational_fitc_loc ?predictive mean variance (*! ETA *)
    end

    module Cov_sampler = struct
      include Common_cov_sampler
      let calc ?predictive means covariances = calc ~loc:variational_fitc_loc ?predictive means covariances (*! ETA *)
    end
  end

  module Variational_FIC = struct
    include Common
    module Model = Variational_model
    module Covariances = FIC_covariances

    module Sampler = struct
      include Common_sampler
      let calc ?predictive mean variance = calc ~loc:variational_fic_loc ?predictive mean variance (*! ETA *)
    end

    module Cov_sampler = struct
      include Common_cov_sampler
      let calc ?predictive means covariances = calc ~loc:variational_fic_loc ?predictive means covariances (*! ETA *)
    end
  end
end


(* Handling derivatives *)

module type Deriv_sig = functor (Spec : Specs.Deriv) ->
  Sigs.Deriv
    with module Eval.Spec = Spec.Eval
    with module Deriv.Spec = Spec

(* Computations shared by FIC and FITC, and standard and variational
   version for derivatives *)
module Make_common_deriv (Spec : Specs.Deriv) = struct
  (* Eval modules *)

  module Eval_common = Make_common (Spec.Eval)
  module Eval_inducing = Eval_common.Inducing
  module Eval_inputs = Eval_common.Inputs
  module Eval_model = Eval_common.Common_model
  module Eval_trained = Eval_common.Trained

  (* Kind of model *)
  type model_kind = Standard | Variational

  module Deriv_common = struct

    (* Derivative modules *)

    module Spec = Spec

    (* Derivative of inducing points *)
    module Inducing = struct
      type ('D, 'd, 'm) t = { eval : ('D, 'd, 'm) Eval_inducing.t; shared_upper : ('D, 'd, 'm, 'd) Spec.Inducing.upper } (*! ITP *)

      let calc kernel eval_inducing =
        let km, shared_upper =
          Spec.Inducing.calc_shared_upper kernel eval_inducing
        in
        {
          eval = Eval_inducing.calc_internal kernel eval_inducing km;
          shared_upper;
        }

      let calc_eval inducing = inducing.eval
      let get_kernel inducing = Eval_inducing.get_kernel inducing.eval
    end

    (* Derivative of inputs *)
    module Inputs = struct
      type ('D, 'd, 'm, 'n) t = { (*! ITP *)
        inducing : ('D, 'd, 'm) Inducing.t; (*! ITP *)
        eval : ('D, 'd, 'm, 'n) Eval_inputs.t; (*! ITP *)
        shared_cross : ('D, 'd, 'm, 'n) Spec.Inputs.cross; (*! ITP *)
      }

      let calc inducing points =
        let kernel = Inducing.get_kernel inducing in
        let knm, shared_cross =
          Spec.Inputs.calc_shared_cross kernel
            ~inputs:points ~inducing:inducing.Inducing.eval.Eval_inducing.points
        in
        let eval =
          Eval_inputs.calc_internal points inducing.Inducing.eval knm
        in
        { inducing; eval; shared_cross }

      let calc_eval t = t.eval

      let get_kernel inputs = Inducing.get_kernel inputs.inducing
    end

    (* Derivative of hyper parameters *)
    module Shared = struct
      type ('D, 'd, 'm, 'n) shared = { (*! ITP *)
        km : ('m, 'm, Slap.cnt) mat; (*! ITP *)
        knm : ('n, 'm, Slap.cnt) mat; (*! ITP *)
        kn_diag : ('n, Slap.cnt) vec; (*! ITP *)
        shared_upper : ('D, 'd, 'm, 'd) Spec.Inducing.upper; (*! ITP *)
        shared_cross : ('D, 'd, 'm, 'n) Spec.Inputs.cross; (*! ITP *)
        shared_diag : ('D, 'd, 'm, 'n) Spec.Inputs.diag; (*! ITP *)
      }

      type ('m, 'n) dfacts = { v_vec : ('n, Slap.cnt) vec; w_mat : ('m, 'm, Slap.cnt) mat; x_mat : ('n, 'm, Slap.cnt) mat } (*! ITP *)
      type ('D, 'd, 'm, 'n) hyper_t = { shared : ('D, 'd, 'm, 'n) shared; dfacts : ('m, 'n) dfacts } (*! ITP *)

      let calc_us_mat eval_model =
        let u_mat = lacpy (Eval_model.get_v_mat eval_model) in
        trsm ~side:Slap.Common.right (*! RF *) ~transa:Slap.Common.trans (*! RF *) u_mat ~a:(Eval_model.get_chol_km eval_model);
        let n = Mat.dim1 u_mat in
        let q_mat = Eval_model.get_q_mat eval_model in
        let s_mat = lacpy (Slap.Mat.submat_dyn n (Mat.dim2 q_mat) q_mat) (*! SUB *) in
        trsm ~side:Slap.Common.right (*! RF *) ~transa:Slap.Common.trans (*! RF *) s_mat ~a:eval_model.Eval_model.r_mat;
        Mat.scal_rows (Eval_model.get_sqrt_is_vec eval_model) s_mat;
        u_mat, s_mat

      let update_tmp tmp v = tmp.x <- tmp.x +. v

      let calc_dkn_diag_term ~v_vec ~kn_diag = function
        | `Vec dkn_diag -> dot ~x:v_vec dkn_diag
        | `Sparse_vec (svec, rows) ->
            check_sparse_vec_sane ~real_n:(Slap.Size.to_int (Vec.dim v_vec)) (*! S2I *) ~svec ~rows;
            let tmp = { x = 0. } in
            for i = 1 to Slap.Size.to_int (Vec.dim svec) (*! S2I *) do
              update_tmp tmp ((Vec.get_dyn v_vec (Slap.Vec.get_dyn rows i)) *. (Vec.get_dyn svec i)) (*! IDX *)
            done;
            tmp.x
        | `Const 0. | `Factor 0. -> 0.
        | `Const c -> c *. Vec.sum v_vec
        | `Factor c -> c *. dot ~x:kn_diag v_vec

      let calc_dkm_term ~w_mat ~km = function
        | `Dense dkm -> Mat.symm2_trace w_mat dkm
        | `Sparse_rows (smat, rows) -> symm2_sparse_trace ~mat:w_mat ~smat ~rows
        | `Const 0. | `Factor 0. | `Diag_const 0. -> 0.
        | `Const c -> c *. sum_symm_mat w_mat
        | `Factor c -> c *. Mat.symm2_trace w_mat km
        | `Diag_vec ddkm ->
            let tmp = { x = 0. } in
            for i = 1 to Slap.Size.to_int (Mat.dim1 w_mat) (*! S2I *) do
              update_tmp tmp ((Vec.get_dyn ddkm i) *. (Mat.get_dyn w_mat i i)) (*! IDX *)
            done;
            tmp.x
        | `Diag_const c ->
            let tmp = { x = 0. } in
            for i = 1 to Slap.Size.to_int (Mat.dim1 w_mat) (*! S2I *) do
              update_tmp tmp (c *. (Mat.get_dyn w_mat i i)) (*! IDX *)
            done;
            tmp.x

      let calc_dknm_term ~x_mat ~knm = function
        | `Dense dknm -> Mat.gemm_trace ~transa:Slap.Common.trans (*! RF *) x_mat ~transb:Slap.Common.normal (*! IF *) dknm
        | `Sparse_cols (sdknm, cols) ->
            let real_n = Slap.Size.to_int (Mat.dim2 x_mat) (*! S2I *) in
            check_sparse_col_mat_sane ~real_n ~smat:sdknm ~cols;
            let m = Slap.Size.to_int (Mat.dim1 sdknm) (*! S2I *) in
            let tmp = { x = 0. } in
            for c = 1 to Slap.Size.to_int (Int_vec.dim cols) (*! S2I *) do
              let real_c = Slap.Vec.get_dyn cols c (*! IDX *) in
              for r = 1 to m do
                update_tmp tmp ((Mat.get_dyn x_mat r real_c) *. (Mat.get_dyn sdknm r c)) (*! IDX *)
              done
            done;
            tmp.x
        | `Const 0. | `Factor 0. -> 0.
        | `Const c -> c *. sum_mat x_mat
        | `Factor c -> c *. Mat.gemm_trace ~transa:Slap.Common.trans (*! RF *) x_mat ~transb:Slap.Common.normal (*! IF *) knm
        | `Sparse_rows (sdknm, rows) ->
            let real_m = Slap.Size.to_int (Mat.dim1 x_mat) (*! S2I *) in
            check_sparse_row_mat_sane ~real_m ~smat:sdknm ~rows;
            let n = Mat.dim2 sdknm in
            let tmp = { x = 0. } in
            for r = 1 to Slap.Size.to_int (Int_vec.dim rows) (*! S2I *) do
              let real_r = Slap.Vec.get_dyn rows r (*! IDX *) in
              for c = 1 to Slap.Size.to_int n (*! S2I *) do
                update_tmp tmp ((Mat.get_dyn x_mat real_r c) *. (Mat.get_dyn sdknm r c)) (*! IDX *)
              done
            done;
            tmp.x

      let calc_log_evidence { shared; dfacts = { v_vec; w_mat; x_mat } } hyper =
        let dkn_diag_term =
          let kn_diag = shared.kn_diag in
          let dkn_diag = Spec.Inputs.calc_deriv_diag shared.shared_diag hyper in
          calc_dkn_diag_term ~v_vec ~kn_diag dkn_diag
        in
        let dkm_term =
          let km = shared.km in
          let dkm = Spec.Inducing.calc_deriv_upper shared.shared_upper hyper in
          calc_dkm_term ~w_mat ~km dkm
        in
        let dknm_term =
          let knm = shared.knm in
          let dknm = Spec.Inputs.calc_deriv_cross shared.shared_cross hyper in
          calc_dknm_term ~x_mat ~knm dknm
        in
        -0.5 *. (dkn_diag_term -. dkm_term) -. dknm_term
    end

    (* Derivative of models *)
    module Common_model = struct

      (* Precomputations for all derivatives *)

      type ('D, 'd, 'm, 'n) t = { (*! ITP *)
        model_kind : model_kind;
        model_shared : ('D, 'd, 'm, 'n) Shared.shared; (*! ITP *)
        eval_model : ('D, 'd, 'm, 'n) Eval_model.t; (*! ITP *)
        inv_km : ('m, 'm, Slap.cnt) mat; (*! ITP *)
        q_diag : ((('n, 'm) Slap.Size.add, 'm) Slap.Size.sub, Slap.cnt) vec; (*! ITP *)
        t_mat : ('m, 'm, Slap.cnt) mat; (*! ITP *)
      }

      let calc_internal model_kind model_shared eval_model inv_km =
        let q_mat = Eval_model.get_q_mat eval_model in
        let n = Slap.Size.sub_dyn (Mat.dim1 q_mat) (Mat.dim2 q_mat) (*! SOP *) in
        let t_mat = lacpy ~uplo:`U inv_km in
        Mat.axpy ~alpha:~-.1. ~x:(ichol eval_model.Eval_model.r_mat) t_mat;
        {
          model_kind; model_shared; eval_model; inv_km; t_mat;
          q_diag = Mat.syrk_diag ~trans:Slap.Common.normal (*! IF *) (Slap.Mat.submat_dyn n (Mat.dim2 q_mat) q_mat) (*! SUB *)
        }

      let calc_common model_kind inputs sigma2 =
        let kernel = Inputs.get_kernel inputs in
        let eval_inputs = inputs.Inputs.eval in
        let kn_diag, shared_diag =
          Spec.Inputs.calc_shared_diag kernel eval_inputs.Eval_inputs.points
        in
        let calc_with_kn_diag =
          match model_kind with
          | Standard -> Eval_model.calc_with_kn_diag
          | Variational -> Eval_common.Variational_model.calc_with_kn_diag
        in
        let eval_model = calc_with_kn_diag eval_inputs sigma2 kn_diag in
        let km = Eval_model.get_km eval_model in
        let knm = Eval_model.get_knm eval_model in
        let kn_diag = Eval_model.get_kn_diag eval_model in
        let chol_km = Eval_model.get_chol_km eval_model in
        let inv_km = ichol chol_km in
        let model_shared =
          {
            Shared.
            km; knm; kn_diag; shared_diag;
            shared_upper = inputs.Inputs.inducing.Inducing.shared_upper;
            shared_cross = inputs.Inputs.shared_cross;
          }
        in
        calc_internal model_kind model_shared eval_model inv_km

      let calc_eval model = model.eval_model
      let calc inputs ~sigma2 = calc_common Standard inputs sigma2

      let update_sigma2 ({ model_kind } as model) sigma2 =
        let update_sigma2 =
          match model_kind with
          | Standard -> Eval_model.update_sigma2
          | Variational -> Eval_common.Variational_model.update_sigma2
        in
        let eval_model = update_sigma2 model.eval_model sigma2 in
        calc_internal model_kind model.model_shared eval_model model.inv_km

      let calc_v1_vec { q_diag; eval_model; model_kind } =
        let is_vec = Eval_model.get_is_vec eval_model in
        let n = Vec.dim is_vec in
        let v1_vec = Vec.create n in
        match model_kind with
        | Standard ->
            for i = 1 to Slap.Size.to_int n (*! S2I *) do Vec.set_dyn v1_vec i ((Vec.get_dyn is_vec i) *. (1. -. (Vec.get_dyn q_diag i))) (*! IDX *) done;
            v1_vec
        | Variational ->
            let r_vec = Eval_model.get_r_vec eval_model in
            for i = 1 to Slap.Size.to_int n (*! S2I *) do
              Vec.set_dyn v1_vec i (*! IDX *)
                ((Vec.get_dyn is_vec i) *. (2. -. (Vec.get_dyn is_vec i) *. (Vec.get_dyn r_vec i) -. (Vec.get_dyn q_diag i))) (*! IDX *)
            done;
            v1_vec

      (* Derivative of sigma2 *)

      let common_calc_log_evidence_sigma2 { eval_model; model_kind } v_vec =
        let sum_v_vec = Vec.sum v_vec in
        let sum =
          match model_kind with
          | Standard -> sum_v_vec
          | Variational -> sum_v_vec -. Vec.sum eval_model.Eval_model.is_vec
        in
        -0.5 *. sum

      let calc_log_evidence_sigma2 model =
        common_calc_log_evidence_sigma2 model (calc_v1_vec model)

      (* Prepare derivative of general hyper-parameters *)

      let prepare_hyper ({ eval_model; t_mat; model_shared } as model) =
        let v_vec = calc_v1_vec model in
        let sqrt_v_vec = Vec.sqrt v_vec in
        let u_mat, x_mat = Shared.calc_us_mat eval_model in
        Mat.scal_rows sqrt_v_vec u_mat;
        let c = lacpy ~uplo:`U t_mat in
        let w_mat = syrk ~trans:Slap.Common.trans (*! RF *) ~alpha:~-.1. u_mat ~beta:1. ~c in
        Mat.scal_rows sqrt_v_vec u_mat;
        Mat.axpy ~alpha:~-.1. ~x:u_mat x_mat;
        let dfacts = { Shared.v_vec; w_mat; x_mat } in
        { Shared.shared = model_shared; dfacts }

      include Shared
    end

    module Cm = Common_model

    module Variational_model = struct
      include Cm

      let calc inputs ~sigma2 = calc_common Variational inputs sigma2
    end

    (* Derivative of trained models *)
    module Trained = struct
      type ('D, 'd, 'm, 'n) t = { (*! ITP *)
        common_model : ('D, 'd, 'm, 'n) Cm.t; (*! ITP *)
        eval_trained : ('D, 'd, 'm, 'n) Eval_trained.t; (*! ITP *)
        w_vec : ('n, Slap.cnt) vec; (*! ITP *)
        v_vec : ('n, Slap.cnt) vec; (*! ITP *)
      }

      let calc common_model ~targets:y =
        let eval_model = common_model.Cm.eval_model in
        let y_, qt_y_ = Eval_trained.prepare_internal eval_model ~y in
        let u_vec = copy y_ in
        let n = Vec.dim y_ in
        let q_mat = Eval_model.get_q_mat eval_model in
        ignore (gemv ~trans:Slap.Common.normal (*! IF *) ~alpha:~-.1. (Slap.Mat.submat_dyn n (Mat.dim2 q_mat) q_mat) qt_y_ ~beta:1. ~y:u_vec); (*! SUB *)
        let l2 = -0.5 *. dot ~x:u_vec y_ in
        let coeffs = qt_y_ in
        trsv ~trans:Slap.Common.normal (*! IF *) eval_model.Eval_model.r_mat coeffs;
        let w_vec = u_vec in
        let sqrt_is_vec = Eval_model.get_sqrt_is_vec eval_model in
        for i = 1 to Slap.Size.to_int n (*! S2I *) do Vec.set_dyn w_vec i ((Vec.get_dyn w_vec i) *. (Vec.get_dyn sqrt_is_vec i)) (*! IDX *) done;
        let v2_vec = Vec.sqr w_vec in
        let v_vec = Cm.calc_v1_vec common_model in
        axpy ~alpha:~-.1. ~x:v2_vec v_vec;
        {
          common_model; w_vec; v_vec;
          eval_trained = Eval_trained.calc_internal eval_model ~y ~coeffs ~l2
        }

      let calc_eval trained = trained.eval_trained

      (* Derivative of sigma2 *)

      let calc_log_evidence_sigma2 { common_model; v_vec } =
        Cm.common_calc_log_evidence_sigma2 common_model v_vec

      (* Derivative of general hyper-parameters *)

      let prepare_hyper { common_model; eval_trained; w_vec; v_vec } =
        let { Cm.eval_model; t_mat; model_shared = shared } = common_model in
        let u_mat, x_mat = Shared.calc_us_mat eval_model in
        let t_vec = eval_trained.Eval_trained.coeffs in
        let w_mat =
          let w_mat = syr ~alpha:~-.1. t_vec (lacpy ~uplo:`U t_mat) in
          let u1_mat = lacpy u_mat in
          Mat.scal_rows (Vec.sqrt (Cm.calc_v1_vec common_model)) u1_mat;
          let w_mat = syrk ~trans:Slap.Common.trans (*! RF *) ~alpha:~-.1. u1_mat ~beta:1. ~c:w_mat in
          let u2_mat = lacpy u_mat ~b:u1_mat in
          Mat.scal_rows w_vec u2_mat;
          syrk ~trans:Slap.Common.trans (*! RF *) u2_mat ~beta:1. ~c:w_mat
        in
        let x_mat =
          Mat.scal_rows v_vec u_mat;
          Mat.axpy ~alpha:~-.1. ~x:u_mat x_mat;
          ger ~alpha:~-.1. w_vec t_vec x_mat
        in
        { Shared.shared; dfacts = { Shared.v_vec; w_mat; x_mat } }

      include Shared
    end

    module Test = struct
      let update_hyper kernel inducing_points points hyper ~eps =
        let value = Spec.Hyper.get_value kernel inducing_points points hyper in
        let value_eps = value +. eps in
        Spec.Hyper.set_values kernel inducing_points points
          [| hyper |] (Vec.make Slap.Size.one (*! SC *) value_eps)

      let is_bad_deriv ~finite_el ~deriv ~tol =
        Float.is_nan finite_el
          || Float.is_nan deriv
          || Float.abs (finite_el -. deriv) > tol

      let check_deriv_hyper ?(eps = 1e-8) ?(tol = 1e-2)
            kernel1 inducing_points1 points1 hyper =
        let kernel2, inducing_points2, points2 =
          update_hyper kernel1 inducing_points1 points1 hyper ~eps
        in
        let eval_inducing1 = Eval_inducing.calc kernel1 inducing_points1 in
        let eval_cross1 = Eval_inputs.calc points1 eval_inducing1 in
        let eval_inducing2 = Eval_inducing.calc kernel2 inducing_points2 in
        let eval_cross2 = Eval_inputs.calc points2 eval_inducing2 in
        let make_finite ~mat1 ~mat2 =
          let res = lacpy mat2 in
          Mat.axpy ~alpha:~-.1. ~x:mat1 res;
          Mat.scal (1. /. eps) res;
          res
        in
        let km1 = eval_inducing1.Eval_inducing.km in
        let finite_dkm =
          make_finite ~mat1:km1 ~mat2:eval_inducing2.Eval_inducing.km
        in
        let inducing1 = Inducing.calc kernel1 inducing_points1 in
        let check_mat ~name ~deriv ~finite ~r ~c =
          let finite_el = Mat.get_dyn finite r c (*! IDX *) in
          if is_bad_deriv ~finite_el ~deriv ~tol then
            failwithf
              "Gpr.Fitc_gp.Make_deriv.Test.check_deriv_hyper: \
              finite difference (%f) and derivative (%f) differ \
              by more than %f on %s.{%d, %d}" finite_el deriv tol name r c ()
        in
        (* Check dkm *)
        begin
          let check = check_mat ~name:"dkm" ~finite:finite_dkm in
          match
            Spec.Inducing.calc_deriv_upper inducing1.Inducing.shared_upper hyper
          with
          | `Dense dkm ->
              let m = Mat.dim1 dkm in
              for c = 1 to Slap.Size.to_int m (*! S2I *) do
                for r = 1 to c do check ~deriv:(Mat.get_dyn dkm r c) (*! IDX *) ~r ~c done
              done
          | `Sparse_rows (sdkm, rows) ->
              let m = Slap.Size.to_int (Int_vec.dim rows) (*! S2I *) in
              let n = Mat.dim2 sdkm in
              let rows_ix_ref = ref 1 in
              for sparse_r = 1 to m do
                let c = Slap.Vec.get_dyn rows sparse_r (*! IDX *) in
                for r = 1 to Slap.Size.to_int n (*! S2I *) do
                  let mat_r, mat_c = if r > c then c, r else r, c in
                  let rows_ix = !rows_ix_ref in
                  let deriv =
                    if
                      rows_ix > m ||
                      let rows_el = Slap.Vec.get_dyn rows rows_ix (*! IDX *) in
                      r < rows_el || c < rows_el
                    then Mat.get_dyn sdkm sparse_r r (*! IDX *)
                    else begin
                      incr rows_ix_ref;
                      Mat.get_dyn sdkm rows_ix c (*! IDX *)
                    end
                  in
                  check ~deriv ~r:mat_r ~c:mat_c
                done;
                rows_ix_ref := 1
              done
          | `Const const ->
              let m = Mat.dim1 km1 in
              for c = 1 to Slap.Size.to_int m (*! S2I *) do
                for r = 1 to c do check ~deriv:const ~r ~c done
              done
          | `Factor const ->
              let m = Mat.dim1 km1 in
              for c = 1 to Slap.Size.to_int m (*! S2I *) do
                for r = 1 to c do check ~deriv:(const *. (Mat.get_dyn km1 r c)) (*! IDX *) ~r ~c done
              done
          | `Diag_vec diag ->
              let m = Mat.dim1 km1 in
              for c = 1 to Slap.Size.to_int m (*! S2I *) do check ~deriv:(Vec.get_dyn diag c) (*! IDX *) ~r:c ~c done
          | `Diag_const const ->
              let m = Mat.dim1 km1 in
              for c = 1 to Slap.Size.to_int m (*! S2I *) do check ~deriv:const ~r:c ~c done
        end;
        (* Check dknm *)
        let inputs = Inputs.calc inducing1 points1 in
        begin
          let knm1 = eval_cross1.Eval_inputs.knm in
          let finite_dknm =
            make_finite ~mat1:knm1 ~mat2:eval_cross2.Eval_inputs.knm
          in
          let check = check_mat ~name:"dknm" ~finite:finite_dknm in
          match
            Spec.Inputs.calc_deriv_cross inputs.Inputs.shared_cross hyper
          with
          | `Dense dknm ->
              let m = Mat.dim1 knm1 in
              for c = 1 to Slap.Size.to_int (Mat.dim2 knm1) (*! S2I *) do
                for r = 1 to Slap.Size.to_int m (*! S2I *) do check ~deriv:(Mat.get_dyn dknm r c) (*! IDX *) ~r ~c done
              done
          | `Sparse_cols (sdknm, cols) ->
              let m = Mat.dim1 sdknm in
              for c = 1 to Slap.Size.to_int (Int_vec.dim cols) (*! S2I *) do
                let real_c = Slap.Vec.get_dyn cols c (*! IDX *) in
                for r = 1 to Slap.Size.to_int m (*! S2I *) do check ~deriv:(Mat.get_dyn sdknm r c) (*! IDX *) ~r ~c:real_c done
              done
          | `Const const ->
              let m = Mat.dim1 knm1 in
              for c = 1 to Slap.Size.to_int (Mat.dim2 knm1) (*! S2I *) do
                for r = 1 to Slap.Size.to_int m (*! S2I *) do check ~deriv:const ~r ~c done
              done
          | `Factor const ->
              let m = Mat.dim1 knm1 in
              for c = 1 to Slap.Size.to_int (Mat.dim2 knm1) (*! S2I *) do
                for r = 1 to Slap.Size.to_int m (*! S2I *) do check ~deriv:(const *. (Mat.get_dyn knm1 r c)) (*! IDX *) ~r ~c done
              done
          | `Sparse_rows (sdknm, rows) ->
              let n = Mat.dim2 sdknm in
              for r = 1 to Slap.Size.to_int (Int_vec.dim rows) (*! S2I *) do
                let real_r = Slap.Vec.get_dyn rows r (*! IDX *) in
                for c = 1 to Slap.Size.to_int n (*! S2I *) do check ~deriv:(Mat.get_dyn sdknm r c) (*! IDX *) ~r:real_r ~c done
              done
        end;
        (* Check dkn diag *)
        begin
          let kn_diag1, shared_diag =
            Spec.Inputs.calc_shared_diag kernel1 points1
          in
          let kn_diag2 = Spec.Eval.Inputs.calc_diag kernel2 points2 in
          let finite_dkn_diag =
            let res = copy kn_diag2 in
            axpy ~alpha:~-.1. ~x:kn_diag1 res;
            scal (1. /. eps) res;
            res
          in
          let check ~deriv ~r =
            let finite_el = Vec.get_dyn finite_dkn_diag r (*! IDX *) in
            if is_bad_deriv ~finite_el ~deriv ~tol then
              failwithf
                "Gpr.Fitc_gp.Make_deriv.Test.check_deriv_hyper: \
                finite difference (%f) and derivative (%f) differ \
                by more than %f on dkn_diag.{%d}" finite_el deriv tol r ()
          in
          match Spec.Inputs.calc_deriv_diag shared_diag hyper with
          | `Vec dkn_diag ->
              for r = 1 to Slap.Size.to_int (Vec.dim dkn_diag) (*! S2I *) do check ~deriv:(Vec.get_dyn dkn_diag r) (*! IDX *) ~r done
          | `Sparse_vec (sdkn_diag, cols) ->
              let n = Int_vec.dim cols in
              for r = 1 to Slap.Size.to_int n (*! S2I *) do check ~deriv:(Vec.get_dyn sdkn_diag r) ~r:(Slap.Vec.get_dyn cols r) (*! IDX *) done
          | `Const const ->
              for r = 1 to Slap.Size.to_int (Vec.dim kn_diag1) (*! S2I *) do check ~deriv:const ~r done
          | `Factor const ->
              for r = 1 to Slap.Size.to_int (Vec.dim kn_diag1) (*! S2I *) do
                check ~deriv:(const *. (Vec.get_dyn kn_diag1 r)) (*! IDX *) ~r
              done
        end

      let self_test ?(eps = 1e-8) ?(tol = 1e-2)
            kernel1 inducing_points1 points1 ~sigma2 ~targets hyper =
        let inducing1 = Inducing.calc kernel1 inducing_points1 in
        let inputs1 = Inputs.calc inducing1 points1 in
        let deriv_model = Cm.calc inputs1 ~sigma2 in
        let eval_model1 = Cm.calc_eval deriv_model in
        let model_log_evidence1 = Eval_model.calc_log_evidence eval_model1 in
        let deriv_trained = Trained.calc deriv_model ~targets in
        let eval_trained1 = Trained.calc_eval deriv_trained in
        let trained_log_evidence1 =
          Eval_trained.calc_log_evidence eval_trained1
        in
        let check ~name ~before ~after ~deriv =
          let finite_el = (after -. before) /. eps in
          if is_bad_deriv ~finite_el ~deriv ~tol then
            failwithf
              "Gpr.Fitc_gp.Make_deriv.Test.self_test: \
              finite difference (%f) and derivative (%f) differ \
              by more than %f on %s" finite_el deriv tol name ()
        in
        match hyper with
        | `Sigma2 ->
            let eval_model2 =
              let sigma2 = sigma2 +. eps in
              Eval_model.calc inputs1.Inputs.eval ~sigma2
            in
            let model_log_evidence2 =
              Eval_model.calc_log_evidence eval_model2
            in
            let model_deriv = Cm.calc_log_evidence_sigma2 deriv_model in
            check ~name:"sigma2(model)"
              ~before:model_log_evidence1 ~after:model_log_evidence2
              ~deriv:model_deriv;
            let eval_trained2 = Eval_trained.calc eval_model2 ~targets in
            let trained_log_evidence2 =
              Eval_trained.calc_log_evidence eval_trained2
            in
            let trained_deriv = Trained.calc_log_evidence_sigma2 deriv_trained in
            check ~name:"sigma2(trained)"
              ~before:trained_log_evidence1 ~after:trained_log_evidence2
              ~deriv:trained_deriv
        | `Hyper hyper ->
            let kernel2, inducing_points2, points2 =
              update_hyper kernel1 inducing_points1 points1 hyper ~eps
            in
            let eval_inducing2 = Eval_inducing.calc kernel2 inducing_points2 in
            let eval_inputs2 = Eval_inputs.calc points2 eval_inducing2 in
            let eval_model2 = Eval_model.calc eval_inputs2 ~sigma2 in
            let model_log_evidence2 =
              Eval_model.calc_log_evidence eval_model2
            in
            let model_hyper_t = Cm.prepare_hyper deriv_model in
            let model_deriv = Cm.calc_log_evidence model_hyper_t hyper in
            check ~name:"hyper(model)"
              ~before:model_log_evidence1 ~after:model_log_evidence2
              ~deriv:model_deriv;
            let eval_trained2 = Eval_trained.calc eval_model2 ~targets in
            let trained_log_evidence2 =
              Eval_trained.calc_log_evidence eval_trained2
            in
            let trained_hyper_t = Trained.prepare_hyper deriv_trained in
            let trained_deriv = Trained.calc_log_evidence trained_hyper_t hyper in
            check ~name:"hyper(trained)"
              ~before:trained_log_evidence1 ~after:trained_log_evidence2
              ~deriv:trained_deriv
    end

    (* Hyper parameter optimization by evidence maximization
       (type II maximum likelihood) *)
    module Optim = struct
      let get_sigma2 targets = function
        | None -> Vec.sqr_nrm2 targets /. float (Slap.Size.to_int (Vec.dim targets)) (*! S2I *)
        | Some sigma2 when sigma2 < 0. ->
            failwithf "Optim.get_sigma2: sigma2 < 0: %f" sigma2 ()
        | Some sigma2 -> sigma2

      module N1000 = Slap.Size.Of_int_dyn(struct let value = 1000 end) (*! O2L[1] *)
      type 'n default_n_rand_inducing = (('n, Slap.Size.ten) Slap.Size.div, N1000.n) Slap.Size.min (*! O2L[1] *)
      let default_n_rand_inducing n_inputs = (*! O2L[1] *)
        Slap.Size.min (Slap.Size.div_dyn n_inputs Slap.Size.ten (*! SC *)) N1000.value (*! SOP,O2L[1] *)

      let get_kernel_inducing ~kernel (*! O2L[2] *) ~n_rand_inducing (*! O2L[1] *) ~inputs = function
        | None ->
            let n_inducing =
              let n_inputs = Spec.Eval.Inputs.get_n_points inputs in
              (* match n_rand_inducing with *) (*! O2L[1] *)
              (* | None -> min (n_inputs / 10) 1000 *) (*! O2L[1] *)
              (* | Some n_rand_inducing -> *) (*! O2L[1] *)
                  if (Slap.Size.to_int n_rand_inducing) (*! S2I *) < 1 then
                    failwithf
                      "Gpr.Fitc_gp.Optim.get_kernel_inducing: \
                      n_rand_inducing (%d) < 1" (Slap.Size.to_int n_rand_inducing) (*! S2I *) ()
                  else if (Slap.Size.to_int n_rand_inducing) > (Slap.Size.to_int n_inputs) (*! S2I *) then
                    failwithf
                      "Gpr.Fitc_gp.Optim.get_kernel_inducing: \
                      n_rand_inducing (%d) > n_inputs (%d)"
                      (Slap.Size.to_int n_rand_inducing) (Slap.Size.to_int n_inputs) (*! S2I *) ()
                  else n_rand_inducing
            in
            (* let kernel = *) (*! O2L[2] *)
            (*   match kernel with *) (*! O2L[2] *)
            (*   | None -> Eval_inputs.create_default_kernel ~n_inducing inputs *) (*! O2L[2] *)
            (*   | Some kernel -> kernel *) (*! O2L[2] *)
            (* in *) (*! O2L[2] *)
            (
              kernel,
              Eval_inducing.choose_n_random_inputs kernel ~n_inducing inputs
            )
        | Some inducing ->
            (* match kernel with *) (*! O2L[2] *)
            (* | None -> *) (*! O2L[2] *)
            (*     let n_inducing = Spec.Eval.Inducing.get_n_points inducing in *) (*! O2L[2] *)
            (*     ( *) (*! O2L[2] *)
            (*       Eval_inputs.create_default_kernel ~n_inducing inputs, *) (*! O2L[2] *)
            (*       inducing *) (*! O2L[2] *)
            (*     ) *) (*! O2L[2] *)
            (* | Some kernel -> *) (*! O2L[2] *) kernel, inducing

      let get_default_hypers ~kernel ~n_rand_inducing ~inputs = (*! EGPT[2] *)
        let kernel, inducing = (*! EGPT[2] *)
          get_kernel_inducing ~kernel (*! O2L[2] *) ~n_rand_inducing ~inputs None (*! EGPT[2] *)
        in (*! EGPT[2] *)
        Spec.Hyper.get_all kernel inducing inputs (*! EGPT[2] *)

      let get_hypers_vals kernel inducing points n_hypers (*! EGPT[2] *) hypers =
        (* let hypers = *) (*! EGPT[2] *)
        (*   match hypers with *) (*! EGPT[2] *)
        (*   | None -> Spec.Hyper.get_all kernel inducing points *) (*! EGPT[2] *)
        (*   | Some hypers -> hypers *) (*! EGPT[2] *)
        (* in *) (*! EGPT[2] *)
        if Array.length hypers <> Slap.Size.to_int n_hypers (*! S2I *) then invalid_arg "get_hypers_vals"; (*! EGPT[2] *)
        (* let n_hypers = Array.length hypers in *)(*! EGPT[2] *)
        let hyper_vals =
          Vec.init n_hypers (fun i1 ->
            Spec.Hyper.get_value kernel inducing points hypers.(i1 - 1))
        in
        hypers, hyper_vals

      module Gsl = struct
        exception Optim_exception of exn

        let check_exception seen_exception_ref res =
          if Float.classify res = Float.Class.Nan then
            match !seen_exception_ref with
            | None ->
                failwith "Gpr.Optim.Gsl: optimization function returned nan"
            | Some exc -> raise (Optim_exception exc)

        let ignore_report ~iter:_ _ = ()

        let train
              ?(step = 1e-1) ?(tol = 1e-1) ?(epsabs = 1e-1)
              ?(report_trained_model = ignore_report)
              ?(report_gradient_norm = ignore_report)
              ~kernel (*! O2L[2] *) ?sigma2 ?inducing ~n_rand_inducing (*! O2L[1] *)
              ?(learn_sigma2 = true) ~n_hypers ~hypers (*! EGPT[2] *) ~inputs ~targets () =
          let sigma2 = get_sigma2 targets sigma2 in
          let kernel, inducing =
            get_kernel_inducing ~kernel (*! O2L[2] *) ~n_rand_inducing (*! O2L[1] *) ~inputs inducing
          in
          let hypers, hyper_vals =
            get_hypers_vals kernel (*! O2L[2] *) inducing inputs n_hypers hypers (*! EGPT[2] *)
          in
          (*let n_hypers = Array.length hypers in*) (*! EGPT[2] *)
          let n_gsl_hypers, gsl_hypers =
            if learn_sigma2 then
              let n_gsl_hypers = 1 + (Slap.Size.to_int n_hypers) (*! S2I *) in
              let gsl_hypers = Gsl.Vector.create n_gsl_hypers in
              gsl_hypers.{0} <- log sigma2;
              for i = 1 to Slap.Size.to_int n_hypers (*! S2I *) do gsl_hypers.{i} <- Vec.get_dyn hyper_vals i (*! IDX *) done;
              n_gsl_hypers, gsl_hypers
            else
              let gsl_hypers = Gsl.Vector.create (Slap.Size.to_int n_hypers) (*! S2I *) in
              for i = 1 to Slap.Size.to_int n_hypers (*! S2I *) do
                gsl_hypers.{i - 1} <- Vec.get_dyn hyper_vals i (*! IDX *)
              done;
              (Slap.Size.to_int n_hypers) (*! S2I *), gsl_hypers
           in
          let module Gd = Gsl.Multimin.Deriv in
          let sigma2_ref = ref sigma2 in
          let update_hypers =
            if learn_sigma2 then
              (fun ~gsl_hypers ->
                sigma2_ref := exp gsl_hypers.{0};
                let hyper_vals = Vec.create n_hypers in
                for i = 1 to Slap.Size.to_int n_hypers (*! S2I *) do Vec.set_dyn hyper_vals i gsl_hypers.{i} (*! IDX *) done;
                Spec.Hyper.set_values kernel inducing inputs hypers hyper_vals)
            else
              (fun ~gsl_hypers ->
                let hyper_vals = Vec.create n_hypers in
                for i = 1 to Slap.Size.to_int n_hypers (*! S2I *) do
                  Vec.set_dyn hyper_vals i gsl_hypers.{i - 1} (*! IDX *)
                done;
                Spec.Hyper.set_values kernel inducing inputs hypers hyper_vals)
          in
          let seen_exception_ref = ref None in
          let wrap_seen_exception f =
            try f () with exc -> seen_exception_ref := Some exc; raise exc
          in
          let best_model_ref = ref None in
          let get_best_model () =
            match !best_model_ref with
            | None -> assert false  (* impossible *)
            | Some (trained, _) -> trained
          in
          let iter_count = ref 1 in
          let update_best_model trained log_evidence =
            match !best_model_ref with
            | Some (_, old_log_evidence)
              when old_log_evidence >= log_evidence -> ()
            | _ ->
                report_trained_model ~iter:!iter_count trained;
                best_model_ref := Some (trained, log_evidence)
          in
          let multim_f ~x:gsl_hypers =
            let kernel, inducing, inputs = update_hypers ~gsl_hypers in
            let eval_inducing = Eval_inducing.calc kernel inducing in
            let eval_inputs = Eval_inputs.calc inputs eval_inducing in
            let model = Eval_model.calc eval_inputs ~sigma2:!sigma2_ref in
            let trained = Eval_trained.calc model ~targets in
            let log_evidence = Eval_trained.calc_log_evidence trained in
            update_best_model trained log_evidence;
            -. log_evidence
          in
          let multim_f ~x = wrap_seen_exception (fun () -> multim_f ~x) in
          let multim_dcommon ~x:gsl_hypers ~g:gradient =
            let kernel, inducing, inputs = update_hypers ~gsl_hypers in
            let deriv_inducing = Inducing.calc kernel inducing in
            let deriv_inputs = Inputs.calc deriv_inducing inputs in
            let dmodel = Cm.calc ~sigma2:!sigma2_ref deriv_inputs in
            let trained = Trained.calc dmodel ~targets in
            if learn_sigma2 then
              let dlog_evidence_dsigma2 =
                Trained.calc_log_evidence_sigma2 trained
              in
              gradient.{0} <- -. dlog_evidence_dsigma2 *. !sigma2_ref;
              if Slap.Size.to_int n_hypers (*! S2I *) = 0 then ()
              else
                let hyper_t = Trained.prepare_hyper trained in
                for i1 = 1 to Slap.Size.to_int n_hypers (*! S2I *) do
                  gradient.{i1} <-
                    -. Trained.calc_log_evidence hyper_t hypers.(i1 - 1)
                done;
            else begin
              let hyper_t = Trained.prepare_hyper trained in
              for i = 0 to (Slap.Size.to_int n_hypers) - 1 (*! S2I *) do
                gradient.{i} <- -. Trained.calc_log_evidence hyper_t hypers.(i)
              done;
            end;
            trained
          in
          let multim_df ~x ~g = ignore (multim_dcommon ~x ~g) in
          let multim_df ~x ~g =
            wrap_seen_exception (fun () -> multim_df ~x ~g)
          in
          let multim_fdf ~x ~g =
            let deriv_trained = multim_dcommon ~x ~g in
            let trained = Trained.calc_eval deriv_trained in
            let log_evidence = Eval_trained.calc_log_evidence trained in
            update_best_model trained log_evidence;
            -. log_evidence
          in
          let multim_fdf ~x ~g =
            wrap_seen_exception (fun () -> multim_fdf ~x ~g)
          in
          let multim_fun_fdf = { Gsl.Fun.multim_f; multim_df; multim_fdf } in
          let mumin =
            Gd.make Gd.VECTOR_BFGS2 n_gsl_hypers
              multim_fun_fdf ~x:gsl_hypers ~step ~tol
          in
          let gsl_dhypers = Gsl.Vector.create n_gsl_hypers in
          let rec loop () =
            let neg_log_likelihood =
              Gd.minimum ~x:gsl_hypers ~g:gsl_dhypers mumin
            in
            check_exception seen_exception_ref neg_log_likelihood;
            let gnorm = Gsl.Blas.nrm2 gsl_dhypers in
            begin
              try report_gradient_norm ~iter:!iter_count gnorm
              with exc -> raise (Optim_exception exc)
            end;
            if gnorm < epsabs then get_best_model ()
            else begin
              incr iter_count;
              Gd.iterate mumin;
              loop ()
            end
          in
          loop ()
      end

      let calc_gradient ~learn_sigma2 ~sigma2 ~n_hypers (*! EGPT[2] *) ~hypers ~trained =
        (* let n_hypers = Array.length hypers in *) (*! EGPT[2] *)
        let n_all_hypers = (snd learn_sigma2) n_hypers in (*! FT[1] *)
        let gradient = Vec.create n_all_hypers in
        if fst learn_sigma2 (*! FT[1] *) then
          let dlog_evidence_dsigma2 =
            Trained.calc_log_evidence_sigma2 trained
          in
          Vec.set_dyn gradient 1 (dlog_evidence_dsigma2 *. sigma2); (*! IDX *)
          if Slap.Size.to_int n_hypers = 0 (*! S2I *) then ()
          else
            let hyper_t = Trained.prepare_hyper trained in
            for i = 0 to (Slap.Size.to_int n_hypers) - 1 (*! S2I *) do
              Vec.set_dyn gradient (i + 2) (Trained.calc_log_evidence hyper_t hypers.(i)) (*! IDX *)
            done
        else begin
          let hyper_t = Trained.prepare_hyper trained in
          for i1 = 1 to Slap.Size.to_int n_hypers (*! S2I *) do
            Vec.set_dyn gradient i1 (Trained.calc_log_evidence hyper_t hypers.(i1 - 1)) (*! IDX *)
          done
        end;
        gradient

      let make_test step gradient_norm get_trained
            ?(epsabs = 0.1) ?max_iter ?(report = ignore) t =
        let max_iter =
          match max_iter with
          | None -> -1
          | Some max_iter when max_iter < 0 ->
              failwith "Optim.SMD.test: max_iter < 0"
          | Some max_iter -> max_iter
        in
        let rec loop n ~best_le ~best ~t =
          if n = 0 || gradient_norm t < epsabs then best
          else
            let new_t = step t in
            let best_le, best =
              let new_trained = get_trained new_t in
              let new_log_evidence =
                Eval_trained.calc_log_evidence new_trained
              in
              if new_log_evidence <= best_le then best_le, best
              else begin
                report new_t;
                new_log_evidence, new_t
              end
            in
            loop (n - 1) ~best_le ~best ~t:new_t
        in
        let best_le = Eval_trained.calc_log_evidence (get_trained t) in
        loop max_iter ~best_le ~best:t ~t

      type ('D, 'd, 'm, 'n, 'ki, 'ko) learn_sigma2 = (*! FT[1] *)
          bool * ('ki Slap.Size.t -> 'ko Slap.Size.t) (*! FT[1] *)
      let learn_sigma2 = (true, Slap.Size.succ) (*! FT[1] *)
      let not_learn_sigma2 = (false, fun n_hypers -> n_hypers) (*! FT[1] *)

      module SGD = struct
        type ('D, 'd, 'm, 'n, 'ki, 'ko) t = { (*! ITP *)
          learn_sigma2 : ('D, 'd, 'm, 'n, 'ki, 'ko) learn_sigma2; (*! FT[1] *)
          hypers : Spec.Hyper.t array;
          tau : float;
          eta : float;
          step : int;
          hyper_vals : ('ki, Slap.cnt) vec; (*! ITP *)
          trained : ('D, 'd, 'm, 'n) Trained.t; (*! ITP *)
          gradient : ('ko, Slap.cnt) vec; (*! ITP *)
          gradient_norm : float;
        }

        let create
              ?(tau = 100.) ?eta0:(eta = 1e-3) ?(step = 0) ~kernel (*! O2L[2] *) ?sigma2 ?inducing
              ~n_rand_inducing (*! O2L[1] *) ~learn_sigma2 (*! FT[1] *) ~n_hypers ~hypers (*! EGPT[2] *)
              ~inputs ~targets () =
          let loc = "Gpr.Fitc_gp.Optim.SGD.create" in
          let fail_neg0 what v =
            if v <= 0. then failwithf "%s: %s (%f) <= 0" loc what v ()
          in
          fail_neg0 "tau" tau;
          fail_neg0 "eta0" eta;
          if step < 0 then failwithf "%s: step (%d) < 0" loc step ();
          let sigma2 = get_sigma2 targets sigma2 in
          let kernel, inducing =
            get_kernel_inducing ~kernel (*! O2L[2] *) ~n_rand_inducing (*! O2L[1] *) ~inputs inducing
          in
          let hypers, hyper_vals =
            get_hypers_vals kernel inducing inputs n_hypers (*! EGPT[2] *) hypers
          in
          let trained =
            let deriv_inducing = Inducing.calc kernel inducing in
            let deriv_inputs = Inputs.calc deriv_inducing inputs in
            let dmodel = Cm.calc ~sigma2 deriv_inputs in
            Trained.calc dmodel ~targets
          in
          let gradient = calc_gradient ~learn_sigma2 ~sigma2 ~n_hypers (*! EGPT[2] *) ~hypers ~trained in
          let gradient_norm = nrm2 gradient in
          {
            learn_sigma2; hypers; tau; eta; step; hyper_vals;
            trained; gradient; gradient_norm
          }

        let step t =
          let
            {
              learn_sigma2; hypers; tau; eta; step;
              hyper_vals = old_hyper_vals;
              trained = old_trained;
              gradient = old_gradient;
              gradient_norm = _gradient_norm;
            } = t
          in
          let old_sigma2, old_input_points, old_inducing, old_kernel =
            let eval_trained = Trained.calc_eval old_trained in
            let eval_model = Eval_trained.get_model eval_trained in
            let input_points = Eval_model.get_input_points eval_model in
            let sigma2 = Eval_model.get_sigma2 eval_model in
            let inducing = Eval_model.get_inducing_points eval_model in
            let kernel = Eval_model.get_kernel eval_model in
            sigma2, input_points, inducing, kernel
          in
          let sigma2, hyper_ix =
            if fst learn_sigma2 (*! FT[1] *) then
              exp (log old_sigma2 +. eta *. (Vec.get_dyn old_gradient 1) (*! IDX *)), 2
            else old_sigma2, 1
          in
          let hyper_vals = copy old_hyper_vals in
          axpy ~alpha:eta ~x:(Slap.Vec.subvec_dyn (Vec.dim hyper_vals) ~ofsx:hyper_ix old_gradient) hyper_vals; (*! SUB *)
          let trained =
            let kernel, inducing, input_points =
              Spec.Hyper.set_values
                old_kernel old_inducing old_input_points hypers hyper_vals
            in
            let deriv_inducing = Inducing.calc kernel inducing in
            let deriv_inputs = Inputs.calc deriv_inducing input_points in
            let dmodel = Cm.calc ~sigma2 deriv_inputs in
            let eval_trained = Trained.calc_eval old_trained in
            let targets = Eval_trained.get_targets eval_trained in
            Trained.calc dmodel ~targets
          in
          let gradient = calc_gradient ~learn_sigma2 ~sigma2 ~n_hypers:(Vec.dim hyper_vals) (*! EGPT[2] *) ~hypers ~trained in
          let gradient_norm = nrm2 gradient in
          {
            t with
            hyper_vals; trained; gradient; gradient_norm;
            eta = tau /. (tau +. float step) *. eta;
            step = step + 1;
          }

        let gradient_norm t = t.gradient_norm
        let get_trained t = Trained.calc_eval t.trained
        let get_eta t = t.eta
        let get_step t = t.step

        let test ?epsabs ?max_iter ?report t = make_test step gradient_norm get_trained ?epsabs ?max_iter ?report t (*! ETA *)
      end

      module SMD = struct
        type ('D, 'd, 'm, 'n, 'ki, 'ko) t = { (*! ITP *)
          learn_sigma2 : ('D, 'd, 'm, 'n, 'ki, 'ko) learn_sigma2; (*! FT[1] *)
          hypers : Spec.Hyper.t array;
          eps : float;
          lambda : float;
          mu : float;
          eta : ('ko, Slap.cnt) vec; (*! ITP *)
          nu : ('ko, Slap.cnt) vec; (*! ITP *)
          hyper_vals : ('ki, Slap.cnt) vec; (*! ITP *)
          trained : ('D, 'd, 'm, 'n) Trained.t; (*! ITP *)
          gradient : ('ko, Slap.cnt) vec; (*! ITP *)
          gradient_norm : float;
        }

        let create
              ?(eps = 1e-8) ?lambda ?mu ?eta0 ?nu0 ~kernel (*! O2L[2] *) ?sigma2 ?inducing
              ~n_rand_inducing (*! O2L[1] *) ~learn_sigma2 (*! FT[1] *) ~n_hypers ~hypers (*! EGPT[2] *)
              ~inputs ~targets () =
          let loc = "Gpr.Fitc_gp.Optim.SMD.create" in
          let lambda =
            match lambda with
            | None -> 0.1
            | Some lambda when lambda < 0. || lambda > 1. ->
                failwithf "%s: violating 0 <= lambda(%f) <= 1" loc lambda ()
            | Some lambda -> lambda
          in
          let mu =
            match mu with
            | None -> 1e-3
            | Some mu when mu < 0. ->
                failwithf "%s: violating 0 <= mu(%f)" loc mu ()
            | Some mu -> mu
          in
          let sigma2 = get_sigma2 targets sigma2 in
          let kernel, inducing =
            get_kernel_inducing ~kernel (*! O2L[2] *) ~n_rand_inducing (*! O2L[1] *) ~inputs inducing
          in
          let hypers, hyper_vals =
            get_hypers_vals kernel inducing inputs n_hypers (*! EGPT[2] *) hypers
          in
          let n_all_hypers =
            (* let n_hypers = Array.length hypers in *) (*! FT[1] *)
            (snd learn_sigma2) n_hypers (*! FT[1] *)
          in
          let eta =
            match eta0 with
            | None -> Vec.make n_all_hypers 1e-3
            | Some eta0 ->
               (* let dim_eta0 = Vec.dim eta0 in *) (*! RMDC *)
               (* if dim_eta0 <> n_all_hypers then *) (*! RMDC *)
               (*   failwithf "%s: dim(eta0) = %d <> n_all_hypers(%d)" *) (*! RMDC *)
               (*     loc dim_eta0 n_all_hypers () *) (*! RMDC *)
               (* else *) (*! RMDC *) begin
                  for i = 1 to Slap.Size.to_int n_all_hypers (*! S2I *) do
                    let eta0_i = Vec.get_dyn eta0 i (*! IDX *) in
                    if eta0_i <= 0. then
                      failwithf "%s: eta0.{%d} < 0: %f" loc i eta0_i ()
                  done;
                  eta0
                end
          in
          let nu =
            match nu0 with
            | None -> Vec.make n_all_hypers 1e-3
            | Some nu0 ->
               (* let dim_nu0 = Vec.dim nu0 in *) (*! RMDC *)
               (* if dim_nu0 <> n_all_hypers then *) (*! RMDC *)
               (*   failwithf "%s: dim(nu0) = %d <> n_all_hypers(%d)" *) (*! RMDC *)
               (*     loc dim_nu0 n_all_hypers () *) (*! RMDC *)
               (* else *) (*! RMDC *) nu0
          in
          let trained =
            let deriv_inducing = Inducing.calc kernel inducing in
            let deriv_inputs = Inputs.calc deriv_inducing inputs in
            let dmodel = Cm.calc ~sigma2 deriv_inputs in
            Trained.calc dmodel ~targets
          in
          let gradient = calc_gradient ~learn_sigma2 ~sigma2 ~n_hypers (*! EGPT[2] *) ~hypers ~trained in
          let gradient_norm = nrm2 gradient in
          {
            learn_sigma2; hypers; eps; lambda; mu; eta; nu;
            hyper_vals; trained; gradient; gradient_norm
          }

        let step t =
          let
            {
              learn_sigma2; hypers; eps; lambda; mu;
              eta = old_eta; nu = old_nu;
              hyper_vals = old_hyper_vals;
              trained = old_trained;
              gradient = old_gradient;
              gradient_norm = _gradient_norm;
            } = t
          in
          let eval_trained = Trained.calc_eval old_trained in
          let targets = Eval_trained.get_targets eval_trained in
          let eval_model = Eval_trained.get_model eval_trained in
          let old_input_points = Eval_model.get_input_points eval_model in
          let old_sigma2 = Eval_model.get_sigma2 eval_model in
          let log_old_sigma2 = log old_sigma2 in
          let old_inducing = Eval_model.get_inducing_points eval_model in
          let old_kernel = Eval_model.get_kernel eval_model in
          let n_hypers = Vec.dim old_hyper_vals in (*! EGPT[2] *)
          let lambda_hessian_nu =
            (* Just an approximation.  Would require algorithmic
               differentiation for practical use if exact Hessian-vector
               product is required. *)
            let calc_grad eps =
              let sigma2, hyper_ofs =
                if fst learn_sigma2 (*! FT[1] *) then
                  exp (log_old_sigma2 +. eps *. (Vec.get_dyn old_nu 1) (*! IDX *)), 1
                else old_sigma2, 0
              in
              let kernel, inducing, input_points =
                let hyper_vals = Vec.create n_hypers in
                for i1 = 1 to Slap.Size.to_int n_hypers (*! S2I *) do
                  Vec.set_dyn hyper_vals i1 (*! IDX *)
                    ((Vec.get_dyn old_hyper_vals i1) +. eps *. (Vec.get_dyn old_nu (i1 + hyper_ofs))) (*! IDX *)
                done;
                Spec.Hyper.set_values
                  old_kernel old_inducing old_input_points hypers hyper_vals
              in
              let deriv_inducing = Inducing.calc kernel inducing in
              let deriv_inputs = Inputs.calc deriv_inducing input_points in
              let dmodel = Cm.calc ~sigma2 deriv_inputs in
              let trained = Trained.calc dmodel ~targets in
              calc_gradient ~learn_sigma2 ~sigma2 ~n_hypers (*! EGPT[2] *) ~hypers ~trained
            in
            let res = Vec.sub (calc_grad eps) (calc_grad (-. eps)) in
            scal (lambda /. (2. *. eps)) res;
            res
          in
          let n_all_hypers = Vec.dim old_gradient in
          let eta = Vec.create n_all_hypers in
          for i = 1 to Slap.Size.to_int n_all_hypers (*! S2I *) do
            Vec.set_dyn eta i (*! IDX *)
              ((Vec.get_dyn old_eta i) *. (*! IDX *)
                max 0.5 (1. +. mu *. (Vec.get_dyn old_gradient i) *. (Vec.get_dyn old_nu i))) (*! IDX *)
          done;
          let sigma2, hyper_ix =
            if fst learn_sigma2 (*! FT[1] *) then
              exp (log_old_sigma2 +. (Vec.get_dyn eta 1) *. (Vec.get_dyn old_gradient 1) (*! IDX *)), 2
            else old_sigma2, 1
          in
          let hyper_vals =
            Vec.add old_hyper_vals
              (Vec.mul (Slap.Vec.subvec_dyn n_hypers eta) (Slap.Vec.subvec_dyn n_hypers ~ofsx:hyper_ix old_gradient)) (*! SUB *)
          in
          let nu =
            Vec.mul old_eta (Vec.add old_gradient lambda_hessian_nu)
          in
          axpy ~alpha:lambda ~x:old_nu nu;
          let trained =
            let kernel, inducing, input_points =
              Spec.Hyper.set_values
                old_kernel old_inducing old_input_points hypers hyper_vals
            in
            let deriv_inducing = Inducing.calc kernel inducing in
            let deriv_inputs = Inputs.calc deriv_inducing input_points in
            let dmodel = Cm.calc ~sigma2 deriv_inputs in
            let eval_trained = Trained.calc_eval old_trained in
            let targets = Eval_trained.get_targets eval_trained in
            Trained.calc dmodel ~targets
          in
          let gradient = calc_gradient ~learn_sigma2 (*! FT[1] *) ~sigma2 ~n_hypers (*! EGPT[2] *) ~hypers ~trained in
          let gradient_norm = nrm2 gradient in
          { t with eta; nu; hyper_vals; trained; gradient; gradient_norm }

        let gradient_norm t = t.gradient_norm
        let get_trained t = Trained.calc_eval t.trained
        let get_eta t = t.eta
        let get_nu t = t.nu

        let test ?epsabs ?max_iter ?report t = make_test step gradient_norm get_trained?epsabs ?max_iter ?report t (*! ETA *)
      end

    end

(*
    module Online = struct
      module SGD = struct
        type foo =
          | Under_capacity of Trained.t
          | Max_capacity of fdsa

        type t = {
          kernel : Spec.Eval.Kernel.t;
          eta : float;
          tau : float;
          foo : foo;
        }
      end

      module SMD = struct
        type t
      end

      type t =
        | Untrained of Spec.Eval.Kernel.t
        | SGD of SGD.t
        | SMD of SMD.t

      let sgd
            ?(capacity = 10)
            ?mean_predictor:_ ?(eta0 = 0.1) ?(tau = 0.) kernel =
        (* TODO: parameter checks *)
        ignore (capacity);
        {
          SGD.
          kernel;
          eta = eta0;
          tau;
          mean_predictor;
        }

      let smd ?(capacity = 10) ?eta0 ?(mu = 0.1) ?(lam = 0.1) kernel =
        (* TODO: parameter checks *)
        ignore (capacity, eta0, mu, lam, kernel);
        (assert false (* XXX *))

      let train_sgd _sgd _input ~target:_ =
        (assert false (* XXX *))

      let train_smd _smd _input ~target:_ =
        (assert false (* XXX *))

      let train t input ~target =
        match t with
        | SGD sgd -> train_sgd sgd input ~target
        | SMD smd -> train_smd smd input ~target

      let calc_mean_predictor = function
        | Untrained kernel ->
            let inputs = Eval.Inputs.create [||] in
            let inducing = Eval.Inputs.create_inducing kernel inputs in
            {
              Eval_common.Mean_predictor.
              inducing;
              coeffs = Vec.empty;
            }
        | SGD _ ->
            (assert false (* XXX *))
        | SMD _ ->
            (assert false (* XXX *))

      let calc_co_variance_predictor _t = (assert false (* XXX *))
    end

*)
  end
end

module Make_FITC_deriv (Spec : Specs.Deriv) = struct
  include Make_common_deriv (Spec)

  module Eval = struct
    include Eval_common

    module Model = Common_model

    module Covariances = FITC_covariances

    module Sampler = struct
      include Common_sampler
      let calc ?predictive mean variance = calc ~loc:fitc_loc ?predictive mean variance (*! ETA *)
    end

    module Cov_sampler = struct
      include Common_cov_sampler
      let calc ?predictive means covariances = calc ~loc:fitc_loc ?predictive means covariances (*! ETA *)
    end
  end

  module Deriv = struct
    include Deriv_common
    module Model = Common_model
  end
end

module Make_FIC_deriv (Spec : Specs.Deriv) = struct
  include Make_common_deriv (Spec)

  module Eval = struct
    include Eval_common

    module Model = Common_model

    module Covariances = FIC_covariances

    module Sampler = struct
      include Common_sampler
      let calc ?predictive mean variance = calc ~loc:fic_loc ?predictive mean variance (*! ETA *)
    end

    module Cov_sampler = struct
      include Common_cov_sampler
      let calc ?predictive means covariances = calc ~loc:fic_loc ?predictive means covariances (*! ETA *)
    end
  end

  module Deriv = struct
    include Deriv_common
    module Model = Common_model
  end
end

module Make_variational_FITC_deriv (Spec : Specs.Deriv) = struct
  include Make_common_deriv (Spec)

  module Eval = struct
    include Eval_common

    module Model = Variational_model

    module Covariances = FITC_covariances

    module Sampler = struct
      include Common_sampler
      let calc ?predictive mean variance = calc ~loc:fitc_loc ?predictive mean variance (*! ETA *)
    end

    module Cov_sampler = struct
      include Common_cov_sampler
      let calc ?predictive means covariances = calc ~loc:fitc_loc ?predictive means covariances (*! ETA *)
    end
  end

  module Deriv = struct
    include Deriv_common
    module Model = Variational_model
  end
end

module Make_variational_FIC_deriv (Spec : Specs.Deriv) = struct
  include Make_common_deriv (Spec)

  module Eval = struct
    include Eval_common

    module Model = Variational_model

    module Covariances = FIC_covariances

    module Sampler = struct
      include Common_sampler
      let calc ?predictive mean variance = calc ~loc:fic_loc ?predictive mean variance (*! ETA *)
    end

    module Cov_sampler = struct
      include Common_cov_sampler
      let calc ?predictive means covariances = calc ~loc:fic_loc ?predictive means covariances (*! ETA *)
    end
  end

  module Deriv = struct
    include Deriv_common
    module Model = Variational_model
  end
end

module Make_deriv (Spec : Specs.Deriv) = struct
  module type Sig = Sigs.Deriv
    with module Eval.Spec = Spec.Eval
    with module Deriv.Spec = Spec

  module Common_deriv = Make_common_deriv (Spec)

  module FITC = struct
    include Common_deriv

    module Eval = struct
      include Eval_common

      module Model = Common_model

      module Covariances = FITC_covariances

      module Sampler = struct
        include Common_sampler
        let calc ?predictive mean variance = calc ~loc:fitc_loc ?predictive mean variance (*! ETA *)
      end

      module Cov_sampler = struct
        include Common_cov_sampler
        let calc ?predictive means covariances = calc ~loc:fitc_loc ?predictive means covariances (*! ETA *)
      end
    end

    module Deriv = struct
      include Deriv_common
      module Model = Common_model
    end
  end

  module FIC = struct
    include Common_deriv

    module Eval = struct
      include Eval_common

      module Model = Variational_model

      module Covariances = FIC_covariances

      module Sampler = struct
        include Common_sampler
        let calc ?predictive mean variance = calc ~loc:fitc_loc ?predictive mean variance (*! ETA *)
      end

      module Cov_sampler = struct
        include Common_cov_sampler
        let calc ?predictive means covariances = calc ~loc:fitc_loc ?predictive means covariances (*! ETA *)
      end
    end

    module Deriv = struct
      include Deriv_common
      module Model = Variational_model
    end
  end

  module Variational_FITC = struct
    include Common_deriv

    module Eval = struct
      include Eval_common

      module Model = Variational_model

      module Covariances = FITC_covariances

      module Sampler = struct
        include Common_sampler
      let calc ?predictive mean variance = calc ~loc:fitc_loc ?predictive mean variance (*! ETA *)
      end

      module Cov_sampler = struct
        include Common_cov_sampler
        let calc ?predictive means covariances = calc ~loc:fitc_loc ?predictive means covariances (*! ETA *)
      end
    end

    module Deriv = struct
      include Deriv_common
      module Model = Variational_model
    end
  end

  module Variational_FIC = struct
    include Common_deriv

    module Eval = struct
      include Eval_common

      module Model = Variational_model

      module Covariances = FIC_covariances

      module Sampler = struct
        include Common_sampler
        let calc ?predictive means covariances = calc ~loc:fic_loc ?predictive means covariances (*! ETA *)
      end

      module Cov_sampler = struct
        include Common_cov_sampler
        let calc ?predictive means covariances = calc ~loc:fic_loc ?predictive means covariances (*! ETA *)
      end
    end

    module Deriv = struct
      include Deriv_common
      module Model = Variational_model
    end
  end
end
