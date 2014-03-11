(* File: block_diag.ml

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

type 'n t = { data : ('n, 'n, Slap.cnt) mat array; n : int } (*! ITP *)

let check_square (i, size) mat =
  (* let m = Mat.dim1 mat in *) (*! RMDC *)
  let n = Mat.dim2 mat in
  i + 1, size + (Slap.Size.to_int n) (*! S2I *)
  (* if m = n then i + 1, size + n *) (*! RMDC *)
  (* else *) (*! RMDC *)
  (*   failwithf *) (*! RMDC *)
  (*     "Block_diag.check_square: matrix at index %d not square: m = %d, n = %d" *) (*! RMDC *)
  (*     i m n () *) (*! RMDC *)

let create data =
  { data; n = snd (Array.fold ~f:check_square ~init:(0, 0) data) }

let copy t = { t with data = Array.map ~f:(fun mat -> lacpy mat) t.data }

let reraise_exc loc i exc =
  Exn.reraisef exc "Block_diag.%s: failed at index %d" loc i  ()

let potrf ?jitter t =
  Array.iteri t.data ~f:(fun i mat ->
    try potrf ?jitter mat with exc -> reraise_exc "potrf" i exc)

let potri ?jitter ?factorize t =
  Array.iteri t.data ~f:(fun i mat ->
    try potri ?jitter ?factorize mat with exc -> reraise_exc "potri" i exc)
