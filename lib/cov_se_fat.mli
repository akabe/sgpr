(* File: cov_se_fat.mli

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

(** {6 Feature-rich ("fat") squared exponential covariance} *)

(** The covariance is defined as:

    [k(x, y) = sf^2 * exp(-1/2 * |(Q_i*P*(x-y))|^2)]

    where [sf^2] is the amplitude, [P] a general [d*D] dimensionality reduction
    matrix ([d << D]), and [Q_i] is a [d*d] diagonal matrix containing all
    multiscales for inducing input number [i].

    Note that multiscales must not get smaller than [0.5] in this framework,
    because the overall length scale is considered to be equal to [1], which
    imposes this mathematical constraint for positive-definiteness.  There is no
    need for a variable global length scale, because the dimensionality
    reduction matrix already generalizes this feature anyway.  Hence an
    unconstrained multiscale parameter [q] is stored as [log(q - 0.5)].

    If [x] and [y] are the same inducing input, then and only then extra noise
    (a different noise level for each inducing input) will be added for
    heteroskedasticity.

    Dimensionality reduction, heteroskedasticity, and multiscales are optional
    features and can be easily turned off by setting the parameters to [None].
*)

open Slap.D (*! RID *)

open Interfaces.Specs

module Params : sig
  type ('D, 'd) tproj (*! FT[2] *)
  val tproj_none : ('d, 'd) tproj (*! FT[2] *)
  val tproj_some : ('D, 'd, Slap.cnt) mat -> ('D, 'd) tproj (*! FT[2] *)

  type ('D, 'd, 'm) params = { (*! ITP *)
    d : 'd Slap.Size.t; (*! ITP *)
    log_sf2 : float;
    tproj : ('D, 'd) tproj; (*! FT[2] *)
    log_hetero_skedasticity : ('m, Slap.cnt) vec option; (*! ITP *)
    log_multiscales_m05 : ('d, 'm, Slap.cnt) mat option; (*! ITP *)
  }

  type ('D, 'd, 'm) t = private ('D, 'd, 'm) params (*! ITP *)

  val create : ('D, 'd, 'm) params -> ('D, 'd, 'm) t (*! ITP *)
end

module Eval :
  Eval
    with type ('D, 'd, 'm) Kernel.params = ('D, 'd, 'm) Params.t (*! ITP *)
    with type ('m, 'n) Inducing.t = ('m, 'n, Slap.cnt) mat (*! ITP *)
    with type 'n Input.t = ('n, Slap.cnt) vec (*! ITP *)
    with type ('m, 'n) Inputs.t = ('m, 'n, Slap.cnt) mat (*! ITP *)


(* Derivatives *)

(* module Proj_hyper : sig type t = private { big_dim : int; small_dim : int } end *)
module Proj_hyper : sig type t = { big_dim : int; small_dim : int } end
module Dim_hyper : sig type t = int end
module Inducing_hyper : sig type t = { ind : int; dim : int } end

module Hyper_repr : sig
  type t =
    [
    | `Log_sf2
    | `Proj of Proj_hyper.t
    | `Log_hetero_skedasticity of Dim_hyper.t
    | `Inducing_hyper of Inducing_hyper.t
    | `Log_multiscale_m05 of Inducing_hyper.t
    ]
end

module Deriv :
  Deriv
    with module Eval = Eval
    with type Hyper.t = Hyper_repr.t
