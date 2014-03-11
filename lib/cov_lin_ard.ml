(* File: cov_lin_ard.ml

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

open Core.Std
open Slap.D (*! RID *)

module Params = struct type ('D, 'd, 'm) t = {id : 'k 'cd . ('D, 'd, 'k, 'cd) Gpr_utils.id; log_ells : ('d, Slap.cnt) vec } (*! ITP,FS[1] *)
                       let create log_ells = {id = (fun x -> x), (fun x -> x); log_ells} (*! FS[1] *) end

module Eval = struct
  module Kernel = struct
    type ('D, 'd, 'm) params = ('D, 'd, 'm) Params.t (*! ITP *)
    type ('D, 'd, 'm) t = { params : ('D, 'd, 'm) params; consts : ('d, Slap.cnt) vec } (*! ITP *)

    let create params =
      let log_ells = params.Params.log_ells in
      let d = Vec.dim log_ells in
      let consts = Vec.create d in
      for i = 1 to Slap.Size.to_int d (*! S2I *) do Vec.set_dyn consts i (exp (-. (Vec.get_dyn log_ells i))) (*! IDX *) done;
      { params; consts }

    let get_params k = k.params
  end

  module Inducing = struct
    type ('m, 'n) t = ('m, 'n, Slap.cnt) mat (*! ITP *)

    let get_n_points = Mat.dim2
    let calc_upper _k inducing = syrk ~trans:Slap.Common.trans (*! RF *) inducing
  end

  module Input = struct
    type 'n t = ('n, Slap.cnt) vec (*! ITP *)

    let calc_ard_input { Kernel.consts } input =
      let d = Vec.dim input in
      let ard_input = Vec.create d in
      for i = 1 to Slap.Size.to_int d (*! S2I *) do Vec.set_dyn ard_input i ((Vec.get_dyn consts i) *. (Vec.get_dyn input i)) (*! IDX *) done;
      ard_input

    let eval k input inducing =
      let input = (fst k.Kernel.params.Params.id) input in (*! FS[1] *)
      gemv ~trans:Slap.Common.trans (*! RF *) inducing (calc_ard_input k input)

    let weighted_eval k input inducing ~coeffs =
      dot ~x:coeffs (eval k input inducing)

    let eval_one { Kernel.consts } input =
      let rec loop res i =
        if i = 0 then res
        else
          let x = (Vec.get_dyn consts i) *. (Vec.get_dyn input i) (*! IDX *) in
          loop (res +. x *. x) (i - 1)
      in
      loop 0. (Slap.Size.to_int (Vec.dim input)) (*! S2I *)
  end

  module Inputs = struct
    type ('m, 'n) t = ('m, 'n, Slap.cnt) mat (*! ITP *)

    let create = Mat.of_col_vecs_dyn (*! RID *)
    let get_n_points = Mat.dim2
    let choose_subset inputs indexes = Gpr_utils.choose_cols inputs indexes

    let calc_ard_inputs { Kernel.params (*! FS[1] *); consts } inputs =
      let inputs = (snd params.Params.id) inputs in (*! FS[1] *)
      let res = lacpy inputs in
      Mat.scal_rows consts res;
      res

    let create_inducing = calc_ard_inputs

    type 'D default_kernel_size = 'D (*! DKS[1] *)

    let create_default_kernel_params inputs ~n_inducing:_ =
      { Params.id = (fun x -> x), (fun x -> x); (*! FS[1] *) log_ells = Vec.make (Mat.dim1 inputs) 0. }

    let calc_upper k inputs = syrk ~trans:Slap.Common.trans (*! RF *) (calc_ard_inputs k inputs)
    let calc_diag k inputs = Mat.syrk_diag ~trans:Slap.Common.trans (*! RF *) (calc_ard_inputs k inputs)

    let calc_cross k ~inputs ~inducing =
      gemm ~transa:Slap.Common.trans (*! RF *) (calc_ard_inputs k inputs) ~transb:Slap.Common.normal (*! IF *) inducing

    let weighted_eval k ~inputs ~inducing ~coeffs =
      gemv ~trans:Slap.Common.normal (*! IF *) (calc_cross k ~inputs ~inducing) coeffs
  end
end

module Deriv = struct
  module Eval = Eval

  module Hyper = struct
    type t = [ `Log_ell of int ]

    let get_all { Eval.Kernel.params } _inducing _inputs =
      Array.init (Slap.Size.to_int (Vec.dim params.Params.log_ells)) (*! S2I *) ~f:(fun d ->
        `Log_ell (d + 1))

    let get_value { Eval.Kernel.params } _inducing _inputs = function
      | `Log_ell i -> Vec.get_dyn params.Params.log_ells i (*! IDX *)

    let set_values k inducing inputs hypers values =
      let { Eval.Kernel.params } = k in
      let log_ells_lazy = lazy (copy params.Params.log_ells) in
      for i = 1 to Array.length hypers do
        match hypers.(i - 1) with
        | `Log_ell d -> Vec.set_dyn (Lazy.force log_ells_lazy) d (Vec.get_dyn values i) (*! IDX *)
      done;
      let new_kernel =
        if Lazy.is_val log_ells_lazy then
          Eval.Kernel.create { Params.id = params.Params.id; (*! FS[1] *) log_ells = Lazy.force log_ells_lazy }
        else k
      in
      new_kernel, inducing, inputs
  end

  module Inducing = struct
    type ('D, 'd, 'm, 'n) upper = ('n, 'm) Eval.Inducing.t (*! ITP,FS[1] *)

    let calc_shared_upper k eval_inducing =
      let upper = Eval.Inducing.calc_upper k eval_inducing in
      upper, eval_inducing

    let calc_deriv_upper _inducing = function
      | `Log_ell _ -> `Const 0.
  end

  module Inputs = struct
    type ('D, 'd, 'm, 'n) diag = ('D, 'd, 'm) Eval.Kernel.t * ('D, 'n) Eval.Inputs.t (*! ITP *)
    type ('D, 'd, 'm, 'n) cross = ('D, 'd, 'm) Eval.Kernel.t * ('D, 'n) Eval.Inputs.t * ('d, 'm) Eval.Inducing.t (*! ITP *)

    let calc_shared_diag k eval_inputs =
      Eval.Inputs.calc_diag k eval_inputs, (k, eval_inputs)

    let calc_shared_cross k ~inputs ~inducing =
      (
        Eval.Inputs.calc_cross k ~inputs ~inducing,
        (k, inputs, inducing)
      )

    let calc_deriv_diag (k, inputs) (`Log_ell d) =
      let n = Mat.dim2 inputs in
      let res = Vec.create n in
      let const = -2. *. (Vec.get_dyn k.Eval.Kernel.consts d) (*! IDX *) in
      for i = 1 to Slap.Size.to_int n (*! S2I *) do
        let el = Mat.get_dyn inputs d i (*! IDX *) in
        Vec.set_dyn res i ((const *. el) *. el) (*! IDX *)
      done;
      `Vec res

    let calc_deriv_cross (k, inputs, inducing) (`Log_ell d) =
      let m = Mat.dim2 inducing in
      let n = Mat.dim2 inputs in
      let res = Mat.create n m in
      let const = -. (Vec.get_dyn k.Eval.Kernel.consts d) (*! IDX *) in
      for c = 1 to Slap.Size.to_int m (*! S2I *) do
        for r = 1 to Slap.Size.to_int n (*! S2I *) do
          Mat.set_dyn res r c (const *. (Mat.get_dyn inducing d c) *. (Mat.get_dyn inputs d r)) (*! IDX *)
        done
      done;
      `Dense res
  end
end
