(* File: cov_se_iso.mli

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

(** {6 Isotropic squared exponential covariance} *)

(** The covariance is defined as:

    [k(x, y) = sf^2 * exp(-1/2 * |1/ell*(x-y)|^2)]

    where [sf^2] is the amplitude, and [ell] is the length scale.
*)

open Slap.D (*! RID *)

open Interfaces.Specs

module Params : sig type ('D, 'd, 'm) t (*! ITP,FS[1] *)
                    val create : log_ell:float -> log_sf2:float -> ('d, 'd, 'm) t (*! FS[1] *) end

type inducing_hyper = { ind : int; dim : int }

module Eval :
  Eval
    with type ('D, 'd, 'm) Kernel.params = ('D, 'd, 'm) Params.t (*! ITP *)
    with type ('m, 'n) Inducing.t = ('m, 'n, Slap.cnt) mat (*! ITP *)
    with type 'n Input.t = ('n, Slap.cnt) vec (*! ITP *)
    with type ('m, 'n) Inputs.t = ('m, 'n, Slap.cnt) mat (*! ITP *)

module Deriv :
  Deriv
    with module Eval = Eval
    with type Hyper.t =
      [ `Log_ell | `Log_sf2 | `Inducing_hyper of inducing_hyper ]
