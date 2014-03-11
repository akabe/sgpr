(* File: cov_se_fat.ml

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

open Interfaces
open Gpr_utils

open Core.Std
open Slap.D (*! RID *)

type ('D, 'd) tnone = { id : 'k 'cd . ('D, 'd, 'k, 'cd) Gpr_utils.id } (*! FT[2] *)
type ('D, 'd) topt = (*! FT[2] *)
  | TNone of ('D, 'd) tnone (*! FT[2] *)
  | TSome of ('D, 'd, Slap.cnt) mat (*! FT[2] *)

let option_map ~f = function None -> None | Some v -> Some (f v)
let option_iter ~f = function None -> () | Some v -> f v
let tproj_iter ~f = function TNone _ -> () | TSome v -> f v (*! FT[2] *)

module Params = struct
  type ('D, 'd) tproj = ('D, 'd) topt (*! FT[2] *)
  let tproj_none = TNone {id = ((fun x -> x), (fun x -> x))} (*! FT[2] *)
  let tproj_some tproj = TSome tproj (*! FT[2] *)

  type ('D, 'd, 'm) params = { (*! ITP *)
    d : 'd Slap.Size.t; (*! ITP *)
    log_sf2 : float;
    tproj : ('D, 'd) tproj; (*! FT[2] *)
    log_hetero_skedasticity : ('m, Slap.cnt) vec option; (*! ITP *)
    log_multiscales_m05 : ('d, 'm, Slap.cnt) mat option; (*! ITP *)
  }

  type ('D, 'd, 'm) t = ('D, 'd, 'm) params (*! ITP *)

  let create (params : ('D, 'd, 'm) params) = (*! ITP *)
    (* let check v_dim name v = *) (*! RMDC *)
    (*   let n = v_dim v in *) (*! RMDC *)
    (*   if n <> params.d then *) (*! RMDC *)
    (*     failwithf *) (*! RMDC *)
    (*       "Cov_se_fat.Params.create: %s projection (%d) disagrees \ *) (*! RMDC *)
    (*       with target dimension d (%d)" name n params.d () *) (*! RMDC *)
    (* in *) (*! RMDC *)
    (* option_iter params.tproj ~f:(check Mat.dim2 "tproj"); *) (*! RMDC *)
    params
end

module Eval = struct
  module Kernel = struct
    type ('D, 'd, 'm) params = ('D, 'd, 'm) Params.t (*! ITP *)

    type ('D, 'd, 'm) t = { (*! ITP *)
      params : ('D, 'd, 'm) params; (*! ITP *)
      sf2 : float;
      hetero_skedasticity : ('m, Slap.cnt) vec option; (*! ITP *)
      multiscales : ('d, 'm, Slap.cnt) mat option; (*! ITP *)
    }

    let create params =
      let hetero_skedasticity =
        option_map params.Params.log_hetero_skedasticity ~f:(Vec.map exp)
      in
      let multiscales =
        let f v = exp v +. 0.5 in
        option_map params.Params.log_multiscales_m05 ~f:(Mat.map f)
      in
      {
        params;
        sf2 = exp params.Params.log_sf2;
        hetero_skedasticity;
        multiscales;
      }

    let get_params k = k.params
  end

  let calc_res_el ~log_sf2 tmp =
    let x = tmp.x in
    tmp.x <- 0.;
    exp (log_sf2 -. 0.5 *. x)

  let calc_upper_vanilla k mat =
    let { Kernel.sf2; params = { Params.d; log_sf2 } } = k in
    let n = Mat.dim2 mat in
    let res = Mat.create n n in
    let tmp = { x = 0. } in
    for c = 1 to Slap.Size.to_int n (*! S2I *) do
      for r = 1 to c - 1 do
        for i = 1 to Slap.Size.to_int d (*! S2I *) do
          let diff = (Mat.get_dyn mat i r) -. (Mat.get_dyn mat i c) (*! IDX *) in
          tmp.x <- tmp.x +. diff *. diff
        done;
        Mat.set_dyn res r c (calc_res_el ~log_sf2 tmp); (*! IDX *)
      done;
      Mat.set_dyn res c c sf2; (*! IDX *)
    done;
    res

  let update_tmp_sum ~tmp ~diff ~scale =
    tmp.x <- tmp.x +. diff *. (diff /. scale) +. log scale

  module Inducing = struct
    type ('m, 'n) t = ('m, 'n, Slap.cnt) mat (*! ITP *)

    let get_n_points = Mat.dim2

    let calc_upper k inducing =
      let m = Mat.dim2 inducing in
      let res =
        match k.Kernel.multiscales with
        | None -> calc_upper_vanilla k inducing
        | Some multiscales ->
            let { Kernel.params = { Params.d; log_sf2 } } = k in
            let res = Mat.create m m in
            let tmp = { x = 0. } in
            for c = 1 to Slap.Size.to_int m (*! S2I *) do
              for r = 1 to c - 1 do
                for i = 1 to Slap.Size.to_int d (*! S2I *) do
                  let diff = (Mat.get_dyn inducing i r) -. (Mat.get_dyn inducing i c) (*! IDX *) in
                  let scale = (Mat.get_dyn multiscales i r) +. (Mat.get_dyn multiscales i c) -. 1. (*! IDX *) in
                  update_tmp_sum ~tmp ~diff ~scale
                done;
                Mat.set_dyn res r c (calc_res_el ~log_sf2 tmp) (*! IDX *)
              done;
              for i = 1 to Slap.Size.to_int d (*! S2I *) do
                let multiscale = Mat.get_dyn multiscales i c (*! IDX *) in
                tmp.x <- tmp.x +. log (multiscale +. multiscale -. 1.)
              done;
              Mat.set_dyn res c c (calc_res_el ~log_sf2 tmp); (*! IDX *)
            done;
            res
      in
      match k.Kernel.hetero_skedasticity with
      | None -> res
      | Some hetero_skedasticity ->
          for i = 1 to Slap.Size.to_int m (*! S2I *) do
            Mat.set_dyn res i i ((Mat.get_dyn res i i) +. (Vec.get_dyn hetero_skedasticity i)) (*! IDX *)
          done;
          res
  end

  module Input = struct
    type 'n t = ('n, Slap.cnt) vec (*! ITP *)

    let eval k input inducing =
      let
        { Kernel.multiscales; params = { Params.d; log_sf2; tproj } } = k
      in
      let projection =
        match tproj with
        | TNone {id = (id, _)} -> id input (*! FT[2] *)
        | TSome tproj -> gemv ~trans:Slap.Common.trans (*! RF *) tproj input (*! FT[2] *)
      in
      let m = Mat.dim2 inducing in
      let res = Vec.create m in
      let tmp = { x = 0. } in
      begin match multiscales with
      | None ->
          for c = 1 to Slap.Size.to_int m (*! S2I *) do
            for i = 1 to Slap.Size.to_int d (*! S2I *) do
              let diff = (Vec.get_dyn projection i) -. (Mat.get_dyn inducing i c) (*! IDX *) in
              tmp.x <- tmp.x +. diff *. diff
            done;
            Vec.set_dyn res c (calc_res_el ~log_sf2 tmp); (*! IDX *)
          done;
      | Some multiscales ->
          for c = 1 to Slap.Size.to_int m (*! S2I *) do
            for i = 1 to Slap.Size.to_int d (*! S2I *) do
              let diff = (Vec.get_dyn projection i) -. (Mat.get_dyn inducing i c) (*! IDX *) in
              let scale = Mat.get_dyn multiscales i c (*! IDX *) in
              update_tmp_sum ~tmp ~diff ~scale
            done;
            Vec.set_dyn res c (calc_res_el ~log_sf2 tmp); (*! IDX *)
          done;
      end;
      res

    let weighted_eval k input inducing ~coeffs =
      dot ~x:(eval k input inducing) coeffs

    let eval_one k _input = k.Kernel.sf2
  end

  module Inputs = struct
    type ('m, 'n) t = ('m, 'n, Slap.cnt) mat (*! ITP *)

    let create = Mat.of_col_vecs_dyn (*! RID *)
    let get_n_points = Mat.dim2
    let choose_subset = choose_cols

    type 'D default_kernel_size = ('D, Slap.Size.ten) Slap.Size.min (*! DKS[1] *)

    let create_default_kernel_params inputs ~n_inducing =
      let big_dim = Mat.dim1 inputs in
      let n_inputs = Mat.dim2 inputs in
      let d = Slap.Size.min big_dim Slap.Size.ten (*! SOP,SC *) in
      let tproj = Mat.create big_dim d in
      let factor = float (Slap.Size.to_int n_inputs) /. float (Slap.Size.to_int big_dim) (*! S2I *) in
      for r = 1 to Slap.Size.to_int big_dim (*! S2I *) do
        let sum_ref = ref 0. in
        for c = 1 to Slap.Size.to_int n_inputs (*! S2I *) do sum_ref := !sum_ref +. (Mat.get_dyn inputs r c) (*! IDX *) done;
        let mean_factor = factor /. !sum_ref in
        for c = 1 to Slap.Size.to_int d (*! S2I *) do
          Mat.set_dyn tproj r c (mean_factor *. (Random.float 2. -. 1.)) (*! IDX *)
        done;
      done;
      {
        Params.
        d;
        log_sf2 = Random.float 2. -. 1.;
        tproj = TSome tproj; (*! FT[2] *)
        log_hetero_skedasticity = Some (Vec.make n_inducing ~-.5.);
        log_multiscales_m05 = Some (Mat.make0 d n_inducing);
      }

    let project k inputs =
      match k.Kernel.params.Params.tproj with
      | TNone {id=(_, id)} -> id inputs (*! FT[2] *)
      | TSome tproj -> gemm ~transa:Slap.Common.trans (*! RF *) tproj ~transb:Slap.Common.normal (*! IF *) inputs (*! FT[2] *)

    let create_inducing = project
    let calc_upper k inputs = calc_upper_vanilla k (project k inputs)
    let calc_diag k inputs = Vec.make (Mat.dim2 inputs) k.Kernel.sf2

    let calc_cross_with_projections k ~projections ~inducing =
      let { Kernel.multiscales; params = { Params.d; log_sf2 } } = k in
      let m = Mat.dim2 inducing in
      let n = Mat.dim2 projections in
      let res = Mat.create n m in
      let tmp = { x = 0. } in
      begin match multiscales with
      | None ->
          for c = 1 to Slap.Size.to_int m (*! S2I *) do
            for r = 1 to Slap.Size.to_int n (*! S2I *) do
              for i = 1 to Slap.Size.to_int d (*! S2I *) do
                let diff = (Mat.get_dyn projections i r) -. (Mat.get_dyn inducing i c) (*! IDX *) in
                tmp.x <- tmp.x +. diff *. diff
              done;
              Mat.set_dyn res r c (calc_res_el ~log_sf2 tmp); (*! IDX *)
            done;
          done;
      | Some multiscales ->
          for c = 1 to Slap.Size.to_int m (*! S2I *) do
            for r = 1 to Slap.Size.to_int n (*! S2I *) do
              for i = 1 to Slap.Size.to_int d (*! S2I *) do
                let diff = (Mat.get_dyn projections i r) -. (Mat.get_dyn inducing i c) (*! IDX *) in
                let scale = Mat.get_dyn multiscales i c (*! IDX *) in
                update_tmp_sum ~tmp ~diff ~scale;
              done;
              Mat.set_dyn res r c (calc_res_el ~log_sf2 tmp); (*! IDX *)
            done;
          done;
      end;
      res

    let calc_cross k ~inputs ~inducing =
      let projections = project k inputs in
      calc_cross_with_projections k ~projections ~inducing

    let weighted_eval k ~inputs ~inducing ~coeffs =
      gemv ~trans:Slap.Common.normal (*! IF *) (calc_cross k ~inputs ~inducing) coeffs
  end
end

module Proj_hyper = struct type t = { big_dim : int; small_dim : int } end
module Dim_hyper = struct type t = int end
module Inducing_hyper = struct type t = { ind : int; dim : int } end

module Hyper_repr = struct
  type t = [
    | `Log_sf2
    | `Proj of Proj_hyper.t
    | `Log_hetero_skedasticity of Dim_hyper.t
    | `Log_multiscale_m05 of Inducing_hyper.t
    | `Inducing_hyper of Inducing_hyper.t
  ]
end

module Deriv = struct
  module Eval = Eval

  module Hyper = struct
    type t = Hyper_repr.t

    let get_all { Eval.Kernel.params } inducing _inputs =
      let
        { Params.d; tproj; log_hetero_skedasticity; log_multiscales_m05 } =
          params
      in
      let m = Mat.dim2 inducing in
      let n_mandatory_hypers = 1 + (Slap.Size.to_int d) * (Slap.Size.to_int m) (*! S2I *) in
      let n_hypers_ref = ref n_mandatory_hypers in
      let update_count_mat maybe_mat =
        option_iter maybe_mat ~f:(fun mat ->
          n_hypers_ref := !n_hypers_ref + Slap.Size.to_int (Mat.dim1 mat) * Slap.Size.to_int (Mat.dim2 mat)) (*! S2I *)
      in
      let update_count_vec maybe_vec =
        option_iter maybe_vec ~f:(fun vec ->
          n_hypers_ref := !n_hypers_ref + Slap.Size.to_int (Vec.dim vec)) (*! S2I *)
      in
      tproj_iter tproj ~f:(fun mat -> (*! FT[2] *)
          n_hypers_ref := !n_hypers_ref + Slap.Size.to_int (Mat.dim1 mat) * Slap.Size.to_int (Mat.dim2 mat)); (*! S2I,FT[2] *)
      update_count_vec log_hetero_skedasticity;
      update_count_mat log_multiscales_m05;
      let n_hypers = !n_hypers_ref in
      let hypers = Array.create ~len:n_hypers `Log_sf2 in
      for ind = 1 to Slap.Size.to_int m (*! S2I *) do
        let indd = (ind - 1) * (Slap.Size.to_int d) (*! S2I *) in
        for dim = 1 to Slap.Size.to_int d (*! S2I *) do
          let inducing_hyper = { Inducing_hyper.ind; dim } in
          hypers.(indd + dim) <- `Inducing_hyper inducing_hyper
        done
      done;
      let pos_ref = ref n_mandatory_hypers in
      tproj_iter (*! FT[2] *) tproj ~f:(fun tproj ->
        let dim = Mat.dim1 tproj in
        for big_dim = 1 to Slap.Size.to_int dim (*! S2I *) do
          for small_dim = 1 to Slap.Size.to_int d (*! S2I *) do
            let pos = !pos_ref in
            pos_ref := pos + 1;
            hypers.(pos) <- `Proj { Proj_hyper.big_dim; small_dim };
          done;
        done);
      option_iter log_hetero_skedasticity ~f:(fun log_hetero_skedasticity ->
        let m = Vec.dim log_hetero_skedasticity in
        for i = 1 to Slap.Size.to_int m (*! S2I *) do
          let pos = !pos_ref in
          pos_ref := pos + 1;
          hypers.(pos) <- `Log_hetero_skedasticity i;
        done);
      option_iter log_multiscales_m05 ~f:(fun log_multiscales_m05 ->
        for ind = 1 to Slap.Size.to_int (Mat.dim2 log_multiscales_m05) (*! S2I *) do
          for dim = 1 to Slap.Size.to_int d (*! S2I *) do
            let pos = !pos_ref in
            pos_ref := pos + 1;
            hypers.(pos) <- `Log_multiscale_m05 { Inducing_hyper.ind; dim };
          done;
        done);
      hypers

    let option_get_value name = function
      | None ->
          failwithf "Deriv.Hyper.option_get_value: %s not supported" name ()
      | Some v -> v
    let tproj_get_value name = function (*! FT[2] *)
      | TNone _ -> failwithf "Deriv.Hyper.topt_get_value: %s not supported" name () (*! FT[2] *)
      | TSome v -> v (*! FT[2] *)

    let get_value { Eval.Kernel.params } inducing _inputs = function
      | `Log_sf2 -> params.Params.log_sf2
      | `Proj { Proj_hyper.big_dim; small_dim } ->
          Mat.get_dyn (tproj_get_value (*! FT[2] *) "tproj" params.Params.tproj) big_dim small_dim (*! IDX *)
      | `Log_hetero_skedasticity dim ->
          Vec.get_dyn (option_get_value "log_hetero_skedasticity" (*! IDX *)
            params.Params.log_hetero_skedasticity) dim (*! IDX *)
      | `Log_multiscale_m05 { Inducing_hyper.ind; dim } ->
          Mat.get_dyn (option_get_value (*! IDX *)
            "log_multiscales_m05" params.Params.log_multiscales_m05) dim ind (*! IDX *)
      | `Inducing_hyper { Inducing_hyper.ind; dim } -> Mat.get_dyn inducing dim ind (*! IDX *)

    let set_values { Eval.Kernel.params } inducing inputs hypers values =
      let log_sf2_ref = ref params.Params.log_sf2 in
      let lazy_opt name f opt_v = lazy (f (option_get_value name opt_v)) in
      let tproj_lazy = lazy (lacpy (tproj_get_value "tproj" params.Params.tproj)) in (*! FT[2] *)
      let log_hetero_skedasticity_lazy =
        lazy_opt "log_hetero_skedasticity"
          copy params.Params.log_hetero_skedasticity
      in
      let log_multiscales_m05_lazy =
        lazy_opt "log_multiscales_m05" lacpy params.Params.log_multiscales_m05
      in
      let inducing_lazy = lazy (lacpy inducing) in
      for i = 1 to Array.length hypers do
        match hypers.(i - 1) with
        | `Log_sf2 -> log_sf2_ref := Vec.get_dyn values i (*! IDX *)
        | `Proj { Proj_hyper.big_dim; small_dim } ->
            Mat.set_dyn (Lazy.force tproj_lazy) big_dim small_dim (Vec.get_dyn values i) (*! IDX *)
        | `Log_hetero_skedasticity dim ->
            Vec.set_dyn (Lazy.force log_hetero_skedasticity_lazy) dim (Vec.get_dyn values i) (*! IDX *)
        | `Log_multiscale_m05 { Inducing_hyper.ind; dim } ->
            Mat.set_dyn (Lazy.force log_multiscales_m05_lazy) dim ind (Vec.get_dyn values i) (*! IDX *)
        | `Inducing_hyper { Inducing_hyper.ind; dim } ->
            Mat.set_dyn (Lazy.force inducing_lazy) dim ind (Vec.get_dyn values i) (*! IDX *)
      done;
      let lift_opt lazy_value value =
        if Lazy.is_val lazy_value then Some (Lazy.force lazy_value)
        else value
      in
      let lift_tproj lazy_value value = (*! FT[2] *)
        if Lazy.is_val lazy_value then TSome (Lazy.force lazy_value) (*! FT[2] *)
        else value in (*! FT[2] *)
      let lift lazy_value value =
        if Lazy.is_val lazy_value then Lazy.force lazy_value
        else value
      in
      let new_kernel =
        Eval.Kernel.create
          {
            Params.
            d = params.Params.d;
            log_sf2 = !log_sf2_ref;
            tproj = (lift_tproj tproj_lazy params.Params.tproj); (*! FT[2] *)
            log_hetero_skedasticity =
              lift_opt log_hetero_skedasticity_lazy
                params.Params.log_hetero_skedasticity;
            log_multiscales_m05 =
              lift_opt
                log_multiscales_m05_lazy params.Params.log_multiscales_m05;
          }
      in
      let new_inducing = lift inducing_lazy inducing in
      new_kernel, new_inducing, inputs
  end

  type ('D, 'd, 'm, 'n) deriv_common = { kernel : ('D, 'd, 'm) Eval.Kernel.t; eval_mat : ('n, 'm, Slap.cnt) mat } (*! ITP *)

  module Inducing = struct
    type ('D, 'd, 'm, 'n) upper = ('n, 'm) Eval.Inducing.t * ('D, 'd, 'm, 'm) deriv_common (*! ITP *)

    let calc_shared_upper kernel inducing =
      let eval_mat = Eval.Inducing.calc_upper kernel inducing in
      eval_mat, (inducing, { kernel; eval_mat })

    let calc_deriv_upper (inducing, { kernel; eval_mat }) hyper =
      match hyper with
      | `Log_sf2 ->
          begin
            match kernel.Eval.Kernel.hetero_skedasticity with
            | None -> `Factor 1.
            | Some hetero_skedasticity ->
                let res = lacpy eval_mat in
                for i = 1 to Slap.Size.to_int (Mat.dim1 res) (*! S2I *) do
                  Mat.set_dyn res i i ((Mat.get_dyn res i i) -. (Vec.get_dyn hetero_skedasticity i)) (*! IDX *)
                done;
                `Dense res
          end
      | `Proj _ -> `Const 0.
      | `Log_hetero_skedasticity dim ->
          begin
            match kernel.Eval.Kernel.hetero_skedasticity with
            | None ->
                failwith (
                    "Cov_se_fat.Deriv.Inducing.calc_deriv_upper: \
                    heteroskedastic modeling disabled, \
                    cannot calculate derivative")
            | Some hetero_skedasticity ->
                let deriv = Vec.make0 (Vec.dim hetero_skedasticity) in
                Vec.set_dyn deriv dim (Vec.get_dyn hetero_skedasticity dim); (*! IDX *)
                (* TODO: sparse diagonal derivatives? *)
                `Diag_vec deriv
          end
      | `Log_multiscale_m05 { Inducing_hyper.ind; dim } ->
          begin match kernel.Eval.Kernel.multiscales with
          | None ->
              failwith (
                  "Cov_se_fat.Deriv.Inducing.calc_deriv_upper: \
                  multiscale modeling disabled, cannot calculate derivative")
          | Some multiscales ->
              let m = Mat.dim2 eval_mat in
              let res = Mat.create Slap.Size.one m (*! SC *) in
              let inducing_dim = Mat.get_dyn inducing dim ind (*! IDX *) in
              let multiscale = Mat.get_dyn multiscales dim ind (*! IDX *) in
              let multiscale_const = multiscale -. 1. in
              let h = 0.5 in
              let multiscale_h = h -. multiscale in
              let multiscale_factor = h *. multiscale_h in
              for i = 1 to ind - 1 do
                let diff = (Mat.get_dyn inducing dim i) -. inducing_dim (*! IDX *) in
                let iscale = 1. /. ((Mat.get_dyn multiscales dim i) +. multiscale_const) (*! IDX *) in
                let sdiff = diff *. iscale in
                let sdiff2 = sdiff *. sdiff in
                let inner = (iscale -. sdiff2) *. multiscale_factor in
                Mat.set_dyn res 1 i (inner *. (Mat.get_dyn eval_mat i ind)) (*! IDX *)
              done;
              begin match kernel.Eval.Kernel.hetero_skedasticity with
              | None ->
                  Mat.set_dyn res 1 ind (*! IDX *)
                    (multiscale_h /. (multiscale +. multiscale_const) (*! IDX *)
                      *. (Mat.get_dyn eval_mat ind ind)); (*! IDX *)
              | Some hetero_skedasticity ->
                  Mat.set_dyn res 1 ind (*! IDX *)
                    (multiscale_h /. (multiscale +. multiscale_const) (*! IDX *)
                      *. ((Mat.get_dyn eval_mat ind ind) -. (Vec.get_dyn hetero_skedasticity ind))); (*! IDX *)
              end;
              for i = ind + 1 to Slap.Size.to_int m (*! S2I *) do
                let diff = (Mat.get_dyn inducing dim i) -. inducing_dim (*! IDX *) in
                let iscale = 1. /. ((Mat.get_dyn multiscales dim i) +. multiscale_const) (*! IDX *) in
                let sdiff = diff *. iscale in
                let sdiff2 = sdiff *. sdiff in
                let inner = (iscale -. sdiff2) *. multiscale_factor in
                Mat.set_dyn res 1 i (inner *. (Mat.get_dyn eval_mat ind i)) (*! IDX *)
              done;
              let rows = Sparse_indices.create Slap.Size.one (*! SC *) in
              Slap.Vec.set_dyn rows 1 ind; (*! IDX *)
              `Sparse_rows (res, rows)
          end
      | `Inducing_hyper { Inducing_hyper.ind; dim } ->
          let m = Mat.dim2 eval_mat in
          let res = Mat.create Slap.Size.one m (*! SC *) in
          let inducing_dim = Mat.get_dyn inducing dim ind (*! IDX *) in
          begin match kernel.Eval.Kernel.multiscales with
          | None ->
              for i = 1 to ind - 1 do
                let diff = (Mat.get_dyn inducing dim i) -. inducing_dim (*! IDX *) in
                Mat.set_dyn res 1 i (diff *. (Mat.get_dyn eval_mat i ind)) (*! IDX *)
              done;
              Mat.set_dyn res 1 ind 0.; (*! IDX *)
              for i = ind + 1 to Slap.Size.to_int m (*! S2I *) do
                let diff = (Mat.get_dyn inducing dim i) -. inducing_dim (*! IDX *) in
                Mat.set_dyn res 1 i (diff *. (Mat.get_dyn eval_mat ind i)) (*! IDX *)
              done
          | Some multiscales ->
              let multiscale_const = (Mat.get_dyn multiscales dim ind) (*! IDX *) -. 1. in
              for i = 1 to ind - 1 do
                let diff = (Mat.get_dyn inducing dim i) (*! IDX *) -. inducing_dim in
                let scale = (Mat.get_dyn multiscales dim i) (*! IDX *) +. multiscale_const in
                Mat.set_dyn res 1 i (diff /. scale *. (Mat.get_dyn eval_mat i ind)) (*! IDX *)
              done;
              Mat.set_dyn res 1 ind 0.; (*! IDX *)
              for i = ind + 1 to Slap.Size.to_int m (*! S2I *) do
                let diff = (Mat.get_dyn inducing dim i (*! IDX *)) -. inducing_dim in
                let scale = (Mat.get_dyn multiscales dim i) (*! IDX *) +. multiscale_const in
                Mat.set_dyn res 1 i (diff /. scale *. (Mat.get_dyn eval_mat ind i)) (*! IDX *)
              done;
          end;
          let rows = Sparse_indices.create Slap.Size.one (*! SC *) in
          Slap.Vec.set_dyn rows 1 ind; (*! IDX *)
          `Sparse_rows (res, rows)
  end

  module Inputs = struct
    (* Diag *)

    type ('D, 'd, 'm, 'n) diag = ('D, 'd, 'm) Eval.Kernel.t (*! ITP,FS[1] *)

    let calc_shared_diag k diag_eval_inputs =
      Eval.Inputs.calc_diag k diag_eval_inputs, k

    let calc_deriv_diag _diag = function
      | `Log_sf2 -> `Factor 1.
      | `Proj _ | `Log_hetero_skedasticity _ | `Log_multiscale_m05 _
      | `Inducing_hyper _ -> `Const 0.

    (* Cross *)

    module Cross = struct
      type ('D, 'd, 'm, 'n) t = { (*! ITP *)
        common : ('D, 'd, 'm, 'n) deriv_common; (*! ITP *)
        inputs : ('D, 'n) Eval.Inputs.t; (*! ITP *)
        inducing : ('d, 'm) Eval.Inducing.t; (*! ITP *)
        projections : ('d, 'n) Eval.Inducing.t; (*! ITP *)
      }
    end

    type ('D, 'd, 'm, 'n) cross = ('D, 'd, 'm, 'n) Cross.t (*! ITP *)

    let calc_shared_cross kernel ~inputs ~inducing =
      let projections = Eval.Inputs.project kernel inputs in
      let eval_mat =
        Eval.Inputs.calc_cross_with_projections kernel ~projections ~inducing
      in
      let shared =
        { Cross.common = { kernel; eval_mat }; inputs; inducing; projections }
      in
      eval_mat, shared

    let check_tproj_available = function
      | TNone _ -> (*! FT[2] *)
          failwith
            "Cov_se_fat.Deriv.Inputs.calc_deriv_cross: \
            tproj disabled, cannot calculate derivative"
      | TSome _ -> () (*! FT[2] *)

    let calc_deriv_cross cross hyper =
      let
        { Cross.common = { kernel; eval_mat }; inputs; inducing; projections } =
          cross
      in
      match hyper with
      | `Log_sf2 -> `Factor 1.
      | `Proj { Proj_hyper.big_dim; small_dim } ->
          check_tproj_available kernel.Eval.Kernel.params.Params.tproj;
          let m = Mat.dim2 inducing in
          let n = Mat.dim2 inputs in
          let res = Mat.create n m in
          begin match kernel.Eval.Kernel.multiscales with
          | None ->
              for c = 1 to Slap.Size.to_int m (*! S2I *) do
                let ind_el = Mat.get_dyn inducing small_dim c (*! IDX *) in
                for r = 1 to Slap.Size.to_int n (*! S2I *) do
                  let alpha = Mat.get_dyn inputs big_dim r (*! IDX *) in
                  let proj = Mat.get_dyn projections small_dim r (*! IDX *) in
                  Mat.set_dyn res r c (alpha *. (ind_el -. proj) *. (Mat.get_dyn eval_mat r c)) (*! IDX *)
                done
              done;
          | Some multiscales ->
              for c = 1 to Slap.Size.to_int m (*! S2I *) do
                let ind_el = Mat.get_dyn inducing small_dim c (*! IDX *) in
                let multiscale = Mat.get_dyn multiscales small_dim c (*! IDX *) in
                for r = 1 to Slap.Size.to_int n (*! S2I *) do
                  let alpha = Mat.get_dyn inputs big_dim r (*! IDX *) in
                  let proj = Mat.get_dyn projections small_dim r (*! IDX *) in
                  Mat.set_dyn res r c (*! IDX *)
                    (alpha *. ((ind_el -. proj) /. multiscale) *. (Mat.get_dyn eval_mat r c)) (*! IDX *)
                done
              done;
          end;
          `Dense res
      | `Log_hetero_skedasticity _ -> `Const 0.
      | `Log_multiscale_m05 { Inducing_hyper.ind; dim } ->
          begin match kernel.Eval.Kernel.multiscales with
          | None ->
              failwith (
                  "Cov_se_fat.Deriv.Inputs.calc_deriv_cross: \
                  multiscale modeling disabled, cannot calculate derivative")
          | Some multiscales ->
            let n = Mat.dim1 eval_mat in
            let res = Mat.create n Slap.Size.one (*! SC *) in
            let inducing_dim = Mat.get_dyn inducing dim ind (*! IDX *) in
            let multiscale = Mat.get_dyn multiscales dim ind (*! IDX *) in
            let h = 0.5 in
            let multiscale_h = h -. multiscale in
            let multiscale_factor = h *. multiscale_h in
            for r = 1 to Slap.Size.to_int n (*! S2I *) do
              let diff = (Mat.get_dyn projections dim r) (*! IDX *) -. inducing_dim in
              let iscale = 1. /. (Mat.get_dyn multiscales dim ind) (*! IDX *) in
              let sdiff = diff *. iscale in
              let sdiff2 = sdiff *. sdiff in
              let inner = (iscale -. sdiff2) *. multiscale_factor in
              Mat.set_dyn res r 1 (inner *. (Mat.get_dyn eval_mat r ind)) (*! IDX *)
            done;
            let cols = Sparse_indices.create Slap.Size.one (*! SC *) in
            Slap.Vec.set_dyn cols 1 ind; (*! IDX *)
            `Sparse_cols (res, cols)
          end
      | `Inducing_hyper { Inducing_hyper.ind; dim } ->
          let n = Mat.dim1 eval_mat in
          let res = Mat.create n Slap.Size.one (*! SC *) in
          let inducing_dim = Mat.get_dyn inducing dim ind (*! IDX *) in
          begin match kernel.Eval.Kernel.multiscales with
          | None ->
              for r = 1 to Slap.Size.to_int n (*! S2I *) do
                let diff = (Mat.get_dyn projections dim r) (*! IDX *) -. inducing_dim in
                Mat.set_dyn res r 1 (diff *. (Mat.get_dyn eval_mat r ind)) (*! IDX *)
              done;
          | Some multiscales ->
              let multiscale_factor = 1. /. (Mat.get_dyn multiscales dim ind) (*! IDX *) in
              for r = 1 to Slap.Size.to_int n (*! S2I *) do
                let diff = (Mat.get_dyn projections dim r) -. inducing_dim (*! IDX *) in
                Mat.set_dyn res r 1 (multiscale_factor *. diff *. (Mat.get_dyn eval_mat r ind)) (*! IDX *)
              done;
          end;
          let cols = Sparse_indices.create Slap.Size.one (*! SC *) in
          Slap.Vec.set_dyn cols 1 ind; (*! IDX *)
          `Sparse_cols (res, cols)
  end
end
