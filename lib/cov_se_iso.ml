(* File: cov_se_iso.ml

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

open Core.Std
open Slap.D (*! RID *)

module Params = struct type ('D, 'd, 'm) t = {id : 'k 'cd . ('D, 'd, 'k, 'cd) Gpr_utils.id; log_ell : float; log_sf2 : float} (*! ITP,FS[1] *)
                       let create ~log_ell ~log_sf2 = {id = (fun x -> x), (fun x -> x); log_ell; log_sf2} (*! FS[1] *) end

type inducing_hyper = { ind : int; dim : int }

module Eval = struct
  module Kernel = struct
    type ('D, 'd, 'm) params = ('D, 'd, 'm) Params.t (*! ITP *)

    type ('D, 'd, 'm) t = { (*! ITP *)
      params : ('D, 'd, 'm) params; (*! ITP *)
      inv_ell2 : float;
      inv_ell2_05 : float;
      log_sf2 : float;
      sf2 : float;
    }

    let create ({ Params.log_sf2; log_ell } as params) =
      let inv_ell2 = exp (-2. *. log_ell) in
      let inv_ell2_05 = -0.5 *. inv_ell2 in
      { params; inv_ell2; inv_ell2_05; log_sf2; sf2 = exp log_sf2 }

    let get_params k = k.params
  end

  open Kernel

  module Inducing = struct
    type ('m, 'n) t = ('m, 'n, Slap.cnt) mat (*! ITP *)

    let get_n_points = Mat.dim2

    let calc_sqr_diff_mat inducing =
      let d = Mat.dim1 inducing in
      let m = Mat.dim2 inducing in
      let res = Mat.create m m in
      let ssqr_diff_ref = ref 0. in
      for c = 1 to Slap.Size.to_int m (*! S2I *) do
        for r = 1 to c - 1 do
          for i = 1 to Slap.Size.to_int d (*! S2I *) do
            let diff = (Mat.get_dyn inducing i c) -. (Mat.get_dyn inducing i r) (*! IDX *) in
            ssqr_diff_ref := !ssqr_diff_ref +. diff *. diff
          done;
          Mat.set_dyn res r c !ssqr_diff_ref; (*! IDX *)
          ssqr_diff_ref := 0.
        done;
        Mat.set_dyn res c c 0. (*! IDX *)
      done;
      res

    let calc_upper_with_sqr_diff_mat k sqr_diff_mat =
      let m = Mat.dim2 sqr_diff_mat in
      let res = Mat.create m m in
      let { inv_ell2_05; log_sf2; sf2 } = k in
      for c = 1 to Slap.Size.to_int m (*! S2I *) do
        for r = 1 to c - 1 do
          Mat.set_dyn res r c (log_sf2 +. inv_ell2_05 *. (Mat.get_dyn sqr_diff_mat r c)); (*! IDX *)
        done;
        Mat.set_dyn res c c sf2; (*! IDX *)
      done;
      res

    let calc_upper k inducing =
      calc_upper_with_sqr_diff_mat k (calc_sqr_diff_mat inducing)
  end

  module Input = struct
    type 'n t = ('n, Slap.cnt) vec (*! ITP *)

    let eval { Kernel.inv_ell2_05; log_sf2 } input inducing =
      let d = Mat.dim1 inducing in
      let m = Mat.dim2 inducing in
      let res = Vec.create m in
      let ssqr_diff_ref = ref 0. in
      for c = 1 to Slap.Size.to_int m (*! S2I *) do
        for i = 1 to Slap.Size.to_int d (*! S2I *) do
          let diff = (Vec.get_dyn input i) -. (Mat.get_dyn inducing i c) (*! IDX *) in
          ssqr_diff_ref := !ssqr_diff_ref +. diff *. diff
        done;
        Vec.set_dyn res c (exp (log_sf2 +. inv_ell2_05 *. !ssqr_diff_ref)); (*! IDX *)
        ssqr_diff_ref := 0.;
      done;
      res

    let weighted_eval k input inducing ~coeffs =
      dot ~x:coeffs (eval k input inducing)

    let eval_one k _input = k.Kernel.sf2
  end

  module Inputs = struct
    type ('m, 'n) t = ('m, 'n, Slap.cnt) mat (*! ITP *)

    let create = Mat.of_col_vecs_dyn (*! RID *)
    let get_n_points = Mat.dim2
    let choose_subset inputs indexes = Gpr_utils.choose_cols inputs indexes
    let create_inducing _kernel inputs = (snd _kernel.Kernel.params.Params.id) inputs (*! FS[1] *)

    type 'D default_kernel_size = 'D (*! DKS[1] *)

    let create_default_kernel_params _inputs ~n_inducing:_ =
      { Params.id = (fun x -> x), (fun x -> x); (*! FS[1] *) log_ell = 0.; log_sf2 = 0. }

    let calc_upper k inputs = Inducing.calc_upper k inputs
    let calc_diag k inputs = Vec.make (Mat.dim2 inputs) k.Kernel.sf2

    let calc_sqr_diff_mat ~inputs ~inducing =
      let d = Mat.dim1 inducing in
      let m = Mat.dim2 inducing in
      let n = Mat.dim2 inputs in
      let res = Mat.create n m in
      let ssqr_diff_ref = ref 0. in
      for c = 1 to Slap.Size.to_int m (*! S2I *) do
        for r = 1 to Slap.Size.to_int n (*! S2I *) do
          for i = 1 to Slap.Size.to_int d (*! S2I *) do
            let diff = (Mat.get_dyn inputs i r) -. (Mat.get_dyn inducing i c) (*! IDX *) in
            ssqr_diff_ref := !ssqr_diff_ref +. diff *. diff
          done;
          Mat.set_dyn res r c !ssqr_diff_ref; (*! IDX *)
          ssqr_diff_ref := 0.
        done
      done;
      res

    let calc_cross_with_sqr_diff_mat k sqr_diff_mat =
      let { Kernel.inv_ell2_05; log_sf2 } = k in
      let n = Mat.dim1 sqr_diff_mat in
      let m = Mat.dim2 sqr_diff_mat in
      let res = Mat.create n m in
      for c = 1 to Slap.Size.to_int m (*! S2I *) do
        for r = 1 to Slap.Size.to_int n (*! S2I *) do
          Mat.set_dyn res r c (exp (log_sf2 +. inv_ell2_05 *. (Mat.get_dyn sqr_diff_mat r c))) (*! IDX *)
        done
      done;
      res

    let calc_cross k ~inputs ~inducing =
      calc_cross_with_sqr_diff_mat k (calc_sqr_diff_mat ~inputs ~inducing)

    let weighted_eval k ~inputs ~inducing ~(coeffs : ('m, _) vec) = (*! ITA[1] *)
      let sqr_diff_mat = calc_sqr_diff_mat ~inputs ~inducing in
      let n = Mat.dim1 sqr_diff_mat in
      let (m : 'm Slap.Size.t) (*! ITA[1] *) = Mat.dim2 sqr_diff_mat in
      (* if Vec.dim coeffs <> m then *) (*! RMDC *)
      (*   failwith "Gpr.Cov_se_iso.Eval.Inputs.weighted_eval: dim(coeffs) <> m"; *) (*! RMDC *)
      let { Kernel.inv_ell2_05; log_sf2 } = k in
      let rec loop r acc c =
        if c = 0 then acc
        else
          let el =
            (Vec.get_dyn coeffs c) *. exp (log_sf2 +. inv_ell2_05 *. (Mat.get_dyn sqr_diff_mat r c)) (*! IDX *)
          in
          loop r (acc +. el) (c - 1)
      in
      Vec.init n (fun r -> loop r 0. (Slap.Size.to_int m)) (*! S2I *)
  end
end

module Deriv = struct
  module Eval = Eval

  type gen_deriv = [ `Log_ell | `Log_sf2 ]

  module Hyper = struct
    type t = [ gen_deriv | `Inducing_hyper of inducing_hyper ]

    let get_all _kernel inducing _inputs =
      let d = Slap.Size.to_int (Mat.dim1 inducing) (*! S2I *) in
      let m = Slap.Size.to_int (Mat.dim2 inducing) (*! S2I *) in
      let n_inducing_hypers = d * m in
      let n_all_hypers = 2 + n_inducing_hypers in
      let hypers = Array.create ~len:n_all_hypers `Log_ell in
      hypers.(1) <- `Log_sf2 ;
      for ind = 1 to m do
        let indd = (ind - 1) * d in
        for dim = 1 to d do
          let inducing_hyper = { ind; dim } in
          hypers.(1 + indd + dim) <- `Inducing_hyper inducing_hyper
        done
      done;
      hypers

    let get_value { Eval.Kernel.params } inducing _inputs = function
      | `Log_ell -> params.Params.log_ell
      | `Log_sf2 -> params.Params.log_sf2
      | `Inducing_hyper { ind; dim } -> Mat.get_dyn inducing dim ind (*! IDX *)

    let set_values { Eval.Kernel.params } inducing inputs hypers values =
      let { Params.log_ell; log_sf2 } = params in
      let log_ell_ref = ref log_ell in
      let log_sf2_ref = ref log_sf2 in
      let inducing_lazy = lazy (lacpy inducing) in
      for i = 1 to Array.length hypers do
        match hypers.(i - 1) with
        | `Log_ell -> log_ell_ref := Vec.get_dyn values i (*! IDX *)
        | `Log_sf2 -> log_sf2_ref := Vec.get_dyn values i (*! IDX *)
        | `Inducing_hyper { ind; dim } ->
            Mat.set_dyn (Lazy.force inducing_lazy) dim ind (Vec.get_dyn values i) (*! IDX *)
      done;
      let new_kernel =
        let log_ell = !log_ell_ref in
        Eval.Kernel.create { Params.id = params.Params.id (*! FS[1] *); log_ell; log_sf2 = !log_sf2_ref }
      in
      let lift lazy_value value =
        if Lazy.is_val lazy_value then Lazy.force lazy_value
        else value
      in
      let new_inducing = lift inducing_lazy inducing in
      new_kernel, new_inducing, inputs
  end

  type ('D, 'd, 'm, 'n) deriv_common = { (*! ITP *)
    kernel : ('D, 'd, 'm) Eval.Kernel.t; (*! ITP *)
    sqr_diff_mat : ('n, 'm, Slap.cnt) mat; (*! ITP *)
    eval_mat : ('n, 'm, Slap.cnt) mat; (*! ITP *)
  }

  module Inducing = struct
    type ('D, 'd, 'm, 'n) upper = ('n, 'm) Eval.Inducing.t * ('D, 'd, 'm, 'm) deriv_common (*! ITP *)

    let calc_shared_upper kernel eval_inducing =
      let module EI = Eval.Inducing in
      let sqr_diff_mat = EI.calc_sqr_diff_mat eval_inducing in
      let eval_mat = EI.calc_upper_with_sqr_diff_mat kernel sqr_diff_mat in
      eval_mat, (eval_inducing, { kernel; sqr_diff_mat; eval_mat })

    let calc_deriv_upper (inducing, common) = function
      | `Log_sf2 -> `Factor 1.
      | `Log_ell ->
          let { sqr_diff_mat; eval_mat; kernel } = common in
          let m = Mat.dim1 sqr_diff_mat in
          let res = Mat.create m m in
          let { Eval.Kernel.inv_ell2 } = kernel in
          for c = 1 to Slap.Size.to_int m (*! S2I *) do
            for r = 1 to c - 1 do
              Mat.set_dyn res r c ((Mat.get_dyn eval_mat r c) *. (Mat.get_dyn sqr_diff_mat r c) *. inv_ell2) (*! IDX *)
            done;
            Mat.set_dyn res c c 0.; (*! IDX *)
          done;
          `Dense res
        | `Inducing_hyper { ind; dim } ->
            let eval_mat = common.eval_mat in
            let m = Mat.dim2 eval_mat in
            let res = Mat.create Slap.Size.one (*! SC *) m in
            let inducing_dim = Mat.get_dyn inducing dim ind (*! IDX *) in
            let inv_ell2 = common.kernel.Eval.Kernel.inv_ell2 in
            for i = 1 to ind - 1 do
              let ind_d = Mat.get_dyn inducing dim i (*! IDX *) in
              Mat.set_dyn res 1 i (*! IDX *)
                (inv_ell2 *. (ind_d -. inducing_dim) *. (Mat.get_dyn eval_mat i ind)) (*! IDX *)
            done;
            Mat.set_dyn res 1 ind 0.; (*! IDX *)
            for i = ind + 1 to Slap.Size.to_int m (*! S2I *) do
              let ind_d = Mat.get_dyn inducing dim i (*! IDX *) in
              Mat.set_dyn res 1 i (*! IDX *)
                (inv_ell2 *. (ind_d -. inducing_dim) *. (Mat.get_dyn eval_mat ind i)) (*! IDX *)
            done;
            let rows = Sparse_indices.create Slap.Size.one (*! SC *) in
            Slap.Vec.set_dyn rows 1 ind; (*! IDX *)
            `Sparse_rows (res, rows)
  end

  module Inputs = struct
    type ('D, 'd, 'm, 'n) diag = ('D, 'd, 'm) Eval.Kernel.t (*! ITP,FS[1] *)

    type ('D, 'd, 'm, 'n) cross = ('D, 'n) Eval.Inputs.t * ('d, 'm) Eval.Inducing.t * ('D, 'd, 'm, 'n) deriv_common (*! ITP *)

    let calc_shared_diag k diag_eval_inputs =
      Eval.Inputs.calc_diag k diag_eval_inputs, k

    let calc_shared_cross kernel ~inputs ~inducing =
      let module EI = Eval.Inputs in
      let sqr_diff_mat = EI.calc_sqr_diff_mat ~inputs ~inducing in
      let eval_mat = EI.calc_cross_with_sqr_diff_mat kernel sqr_diff_mat in
      let shared = inputs, inducing, { kernel; sqr_diff_mat; eval_mat } in
      eval_mat, shared

    let calc_deriv_diag _diag = function
      | `Log_sf2 -> `Factor 1.
      | `Log_ell | `Inducing_hyper _ -> `Const 0.

    let calc_deriv_cross (inputs, inducing, common) = function
      | `Log_sf2 -> `Factor 1.
      | `Log_ell ->
          let { sqr_diff_mat; eval_mat; kernel } = common in
          let n = Mat.dim1 sqr_diff_mat in
          let m = Mat.dim2 sqr_diff_mat in
          let res = Mat.create n m in
          let { Eval.Kernel.inv_ell2 } = kernel in
          for c = 1 to Slap.Size.to_int m (*! S2I *) do
            for r = 1 to Slap.Size.to_int n (*! S2I *) do
              Mat.set_dyn res r c ((Mat.get_dyn eval_mat r c) *. (Mat.get_dyn sqr_diff_mat r c) *. inv_ell2) (*! IDX *)
            done
          done;
          `Dense res
      | `Inducing_hyper { ind; dim } ->
          let eval_mat = common.eval_mat in
          let n = Mat.dim1 eval_mat in
          let res = Mat.create n Slap.Size.one (*! SC *) in
          let indx_d = Mat.get_dyn inducing dim ind (*! IDX *) in
          let inv_ell2 = common.kernel.Eval.Kernel.inv_ell2 in
          for r = 1 to Slap.Size.to_int n (*! S2I *) do
            let inp_d = Mat.get_dyn inputs dim r (*! IDX *) in
            Mat.set_dyn res r 1 (inv_ell2 *. (inp_d -. indx_d) *. (Mat.get_dyn eval_mat r ind)) (*! IDX *)
          done;
          let cols = Sparse_indices.create Slap.Size.one (*! SC *) in
          Slap.Vec.set_dyn cols 1 ind; (*! IDX *)
          `Sparse_cols (res, cols)
  end
end
