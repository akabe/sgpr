(* File: cov_const.ml

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

module Params = struct type ('D, 'd, 'm) t = { log_theta : float } (*! ITP,FS[1] *)
                       let create log_theta = { log_theta } (*! FS[1] *) end

module Eval = struct
  module Kernel = struct
    type ('D, 'd, 'm) params = ('D, 'd, 'm) Params.t (*! ITP *)
    type ('D, 'd, 'm) t = { params : ('D, 'd, 'm) params; const : float } (*! ITP *)

    let create params = { params; const = exp (-2. *. params.Params.log_theta) }

    let get_params k = k.params
  end

  module Inducing = struct
    type ('m, 'n) t = 'n Slap.Size.t (*! ITP,FS[1] *)

    let calc_upper k m = Mat.make m m k.Kernel.const
    let get_n_points m = m
  end

  module Input = struct
    type 'n t = unit (*! ITP,FS[1] *)

    let eval k () m = Vec.make m k.Kernel.const
    let weighted_eval k () _ ~coeffs = k.Kernel.const *. Vec.sum coeffs
    let eval_one k () = k.Kernel.const
  end

  module Inputs = struct
    type ('m, 'n) t = 'n Slap.Size.t (*! ITP,FS[1] *)

    let create _ n a = if Slap.Size.to_int n (*! S2I *) = Array.length a then n else invalid_arg "Gpr.Cov_const.Eval.Inputs.create" (*! EGPT[1] *)
    let get_n_points n = n
    let choose_subset _inputs indexes = Slap.Vec.dim indexes (*! RID *)
    let create_inducing _kernel n = n

    type 'D default_kernel_size = 'D (*! DKS[1] *)

    let create_default_kernel_params _inputs ~n_inducing:_ =
      { Params.log_theta = 0. }

    let calc_upper = Inducing.calc_upper
    let calc_diag k n = Vec.make n k.Kernel.const
    let calc_cross k ~inputs:n ~inducing:m = Mat.make n m k.Kernel.const

    let weighted_eval k ~inputs:_ ~inducing:_ ~coeffs =
      let res = copy coeffs in
      scal k.Kernel.const res;
      res
  end
end

module Deriv = struct
  module Eval = Eval

  module Hyper = struct
    type t = [ `Log_theta ]

    let get_all _kernel _inducing _inputs = [| `Log_theta |]

    let get_value { Eval.Kernel.params } _inducing _inputs = function
      | `Log_theta -> params.Params.log_theta

    let set_values kernel inducing inputs hypers values =
      let { Eval.Kernel.params } = kernel in
      let log_theta_ref = ref params.Params.log_theta in
      let kernel_changed_ref = ref false in
      for i = 1 to Array.length hypers do
        match hypers.(i - 1) with
        | `Log_theta -> log_theta_ref := (Vec.get_dyn values i) (*! IDX *); kernel_changed_ref := true
      done;
      let new_kernel =
        if !kernel_changed_ref then
          Eval.Kernel.create { Params.log_theta = !log_theta_ref }
        else kernel
      in
      new_kernel, inducing, inputs
  end

  let calc_const_deriv k = -2. *. k.Eval.Kernel.const

  module Inducing = struct
    type ('D, 'd, 'm, 'n) upper = { m : int; deriv_const : float } (*! ITP,FS[1] *)

    let calc_shared_upper k m =
      Eval.Inducing.calc_upper k m, { m = Slap.Size.to_int m (*! S2I *); deriv_const = calc_const_deriv k }

    let calc_deriv_upper shared `Log_theta = `Const shared.deriv_const
  end

  module Inputs = struct
    type ('D, 'd, 'm, 'n) diag = { diag_eval_inputs : ('D, 'n) Eval.Inputs.t; diag_const_deriv : float } (*! ITP,FS[1] *)
    type ('D, 'd, 'm, 'n) cross = { cross_const_deriv : float } (*! ITP,FS[1] *)

    let calc_shared_diag k diag_eval_inputs =
      (
        Eval.Inputs.calc_diag k diag_eval_inputs,
        { diag_eval_inputs; diag_const_deriv = calc_const_deriv k }
      )

    let calc_shared_cross k ~inputs ~inducing =
      (
        Eval.Inputs.calc_cross k ~inputs ~inducing,
        { cross_const_deriv = calc_const_deriv k }
      )

    let calc_deriv_diag diag `Log_theta = `Const diag.diag_const_deriv
    let calc_deriv_cross cross `Log_theta = `Const cross.cross_const_deriv
  end
end
