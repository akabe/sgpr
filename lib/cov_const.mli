(* File: cov_const.mli

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

(** {6 Covariance of a constant function} *)

(** The covariance is defined as:

    [k(x, y) = 1/s^2]
    [logtheta = log(s)]
*)

open Interfaces.Specs

module Params : sig type ('D, 'd, 'm) t (*! ITP,FS[1] *)
                    val create : float -> ('d, 'd, 'm) t (*! FS[1] *) end

module Eval :
  Eval
    with type ('D, 'd, 'm) Kernel.params = ('D, 'd, 'm) Params.t (*! ITP *)
    with type ('m, 'n) Inducing.t = 'n Slap.Size.t (*! ITP *)
    with type 'n Input.t = unit (*! ITP *)
    with type ('m, 'n) Inputs.t = 'n Slap.Size.t (*! ITP *)

module Deriv :
  Deriv
    with module Eval = Eval
    with type Hyper.t = [ `Log_theta ]
