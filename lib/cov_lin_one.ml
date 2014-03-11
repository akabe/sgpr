(* File: cov_lin_one.ml

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

module Params = struct type ('D, 'd, 'm) t = {id : 'k 'cd . ('D, 'd, 'k, 'cd) Gpr_utils.id; log_theta : float} (*! ITP,FS[1] *)
                       let create log_theta = {id = (fun x -> x), (fun x -> x); log_theta} (*! FS[1] *) end

module Eval = struct
  module Kernel = struct
    type ('D, 'd, 'm) params = ('D, 'd, 'm) Params.t (*! ITP *)
    type ('D, 'd, 'm) t = { params : ('D, 'd, 'm) params; const : float } (*! ITP *)

    let create params =
      { params; const = exp (-2. *. params.Params.log_theta) }

    let get_params k = k.params
  end

  module Inducing = struct
    type ('m, 'n) t = ('m, 'n, Slap.cnt) mat (*! ITP *)

    let get_n_points = Mat.dim2

    let calc_upper { Kernel.const = alpha } inducing =
      let m = Mat.dim2 inducing in
      syrk ~alpha ~trans:Slap.Common.trans (*! RF *) inducing ~beta:1. ~c:(Mat.make m m alpha)
  end

  module Input = struct
    type 'n t = ('n, Slap.cnt) vec (*! ITP *)

    let eval { Kernel.params (*! FS[1] *); const = alpha } input inducing =
      let input = (fst params.Params.id) input in (*! FS[1] *)
      gemv ~alpha ~trans:Slap.Common.trans (*! RF *) inducing input
        ~beta:1. ~y:(Vec.make (Mat.dim2 inducing) alpha)

    let weighted_eval k input inducing ~coeffs =
      dot ~x:coeffs (eval k input inducing)

    let eval_one k input = k.Kernel.const *. (Vec.sqr_nrm2 input +. 1.)
  end

  module Inputs = struct
    type ('m, 'n) t = ('m, 'n, Slap.cnt) mat (*! ITP *)

    let create = Mat.of_col_vecs_dyn (*! RID *)
    let get_n_points = Mat.dim2
    let choose_subset inputs indexes = Gpr_utils.choose_cols inputs indexes
    let create_inducing _kernel inputs = (snd _kernel.Kernel.params.Params.id) inputs (*! FS[1] *)

    type 'D default_kernel_size = 'D (*! DKS[1] *)

    let create_default_kernel_params _inputs ~n_inducing:_ =
      { Params.id = (fun x -> x), (fun x -> x) (*! FS[1] *); log_theta = 0. }

    let calc_upper = Inducing.calc_upper

    let calc_diag { Kernel.const = alpha } inputs =
      Mat.syrk_diag ~alpha ~trans:Slap.Common.trans (*! RF *) inputs ~beta:1.
        ~y:(Vec.make (Mat.dim2 inputs) alpha)

    let calc_cross { Kernel.params (*! FS[1] *); const = alpha } ~inputs ~inducing =
      let inputs = (snd params.Params.id) inputs in (*! FS[1] *)
      let m = Mat.dim2 inducing in
      let n = Mat.dim2 inputs in
      gemm ~alpha ~transa:Slap.Common.trans (*! RF *) inputs ~transb:Slap.Common.normal (*! IF *) inducing ~beta:1. ~c:(Mat.make n m alpha)

    let weighted_eval k ~inputs ~inducing ~coeffs =
      gemv ~trans:Slap.Common.normal (*! IF *) (calc_cross k ~inputs ~inducing) coeffs
  end
end

module Deriv = struct
  module Eval = Eval

  module Hyper = struct
    type t = [ `Log_theta ]

    let get_all _kernel _inducing _inputs = [| `Log_theta |]

    let get_value { Eval.Kernel.params = params } _inducing _inputs =
      function `Log_theta -> params.Params.log_theta

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
          Eval.Kernel.create { Params.id = params.Params.id (*! FS[1] *); log_theta = !log_theta_ref }
        else kernel
      in
      new_kernel, inducing, inputs
  end

  let calc_deriv_common () `Log_theta = `Factor ~-.2.

  module Inducing = struct
    type ('D, 'd, 'm, 'n) upper = unit (*! ITP,FS[1] *)

    let calc_shared_upper k inducing = Eval.Inducing.calc_upper k inducing, ()
    let calc_deriv_upper = calc_deriv_common
  end

  module Inputs = struct
    type ('D, 'd, 'm, 'n) diag = unit (*! ITP,FS[1] *)
    type ('D, 'd, 'm, 'n) cross = unit (*! ITP,FS[1] *)

    let calc_shared_diag k diag_eval_inputs =
      Eval.Inputs.calc_diag k diag_eval_inputs, ()

    let calc_shared_cross k ~inputs ~inducing =
      Eval.Inputs.calc_cross k ~inputs ~inducing, ()

    let calc_deriv_diag = calc_deriv_common
    let calc_deriv_cross = calc_deriv_common
  end
end
