(* File: utils.ml

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
open Bigarray
open Slap.D (*! RID *)

(* Global definitions *)

module Int_vec = struct
  type ('n, 'cnt_or_dsc) t = ('n, 'cnt_or_dsc) Slap.Common.int_vec (*! ITP *)
  let create n = Slap.Common.create_int_vec n (*! RID *)
  let dim t = Slap.Vec.dim t (*! RID *)
  let sub t n = Slap.Vec.subvec_dyn n t (*! RID *)
end

let debug = ref false
let cholesky_jitter = ref 1e-6

type fast_float_ref = { mutable x : float }

let pi = 4. *. atan 1.
let log_2pi = log (pi +. pi)

let default_rng = Gsl.Rng.make (Gsl.Rng.default ())


(* Testing and I/O functionality *)

let print_int name n = Format.printf "%s: @[%d@]@.@." name n
let print_float name n = Format.printf "%s: @[%.9f@]@.@." name n
let print_vec name vec = Format.printf "%s: @[%a@]@.@." name pp_vec vec
let print_mat name mat = Format.printf "%s: @[%a@]@.@." name pp_mat mat

let timing name f =
  let t1 = Unix.times () in
  let res = f () in
  let t2 = Unix.times () in
  Format.printf "%s %.2f@." name (t2.Unix.tms_utime -. t1.Unix.tms_utime);
  res

(* General matrix functions *)

(* Choose columns of a matrix *)
let choose_cols mat indexes =
  let m = Mat.dim1 mat in
  let n = Slap.Size.to_int (Mat.dim2 mat) (*! S2I *) in
  let k = Int_vec.dim indexes in
  let res = Mat.create m k in
  for c = 1 to Slap.Size.to_int k (*! S2I *) do
    let real_c = Slap.Vec.get_dyn indexes c (*! IDX *) in
    if real_c < 1 || real_c > n then
      failwithf
        "Gpr.Gpr_utils.choose_cols: violating 1 <= index (%d) <= dim (%d)"
        real_c n ()
    else for r = 1 to Slap.Size.to_int m (*! S2I *) do Mat.set_dyn res r c (Mat.get_dyn mat r real_c) (*! IDX *) done
  done;
  res

(* Compute the sum of all elements in a matrix *)
let sum_mat mat = Vec.sum (Mat.as_vec mat)

(* Compute the sum of all elements in a symmetric matrix *)
let sum_symm_mat mat =
  let diag_ref = ref 0. in
  let rest_ref = ref 0. in
  let n = Mat.dim1 mat in
  for c = 1 to Slap.Size.to_int n (*! S2I *) do
    for r = 1 to c - 1 do rest_ref := !rest_ref +. (Mat.get_dyn mat r c) (*! IDX *) done;
    diag_ref := !diag_ref +. (Mat.get_dyn mat c c) (*! IDX *)
  done;
  let rest = !rest_ref in
  rest +. !diag_ref +. rest

(* Computes logarithm of determinant; assumes Cholesky factorized matrix *)
let log_det (mat : ('n, 'n, 'cd) mat) (*! ITA[2] *) =
  let n = Mat.dim1 mat in
  (* if Mat.dim2 mat <> n then failwith "log_det: not a square matrix"; *) (*! RMDC *)
  let rec loop acc i =
    if i = 0 then acc +. acc
    else loop (acc +. log (Mat.get_dyn mat i i) (*! IDX *)) (i - 1)
  in
  loop 0. (Slap.Size.to_int n) (*! S2I *)

(* Solve triangular system *)
let solve_tri ~trans (*! IF *) chol mat =
  let ichol_mat = lacpy mat in
  trtrs ~trans (*! IF *) chol ichol_mat;
  ichol_mat

(* Compute the inverse of a matrix using the cholesky factor *)
let ichol chol =
  let inv = lacpy ~uplo:`U chol in
  potri ~factorize:false inv;
  inv

(* Sparse matrices and vectors *)

(* Checks whether a sparse row matrix is sane *)
let check_sparse_row_mat_sane ~real_m ~smat ~rows =
  if !debug then begin
    if real_m < 0 then
      failwith "Gpr.Gpr_utils.check_sparse_row_mat_sane: real_m < 0";
    let (m : 'n Slap.Size.t) (*! ITA[3] *) = Mat.dim1 smat in
    let (n_rows : 'n Slap.Size.t) (*! ITA[3] *) = Int_vec.dim rows in
    (* if n_rows <> m then *) (*! RMDC *)
    (*   failwithf *) (*! RMDC *)
    (*     "Gpr.Gpr_utils.check_sparse_row_mat_sane: number of rows in \ *) (*! RMDC *)
    (*     sparse matrix (%d) disagrees with size of row array (%d)" *) (*! RMDC *)
    (*     m n_rows (); *) (*! RMDC *)
    let rec loop ~i ~limit =
      if i > 0 then
        let rows_i = Slap.Vec.get_dyn rows i (*! IDX *) in
        if rows_i <= 0 then
          failwithf
            "Gpr.Gpr_utils.check_sparse_row_mat_sane: sparse row %d contains \
            illegal negative real row index %d" i rows_i ()
        else if rows_i > limit then
          failwithf
            "Gpr.Gpr_utils.check_sparse_row_mat_sane: sparse row %d \
            associated with real row index %d violates consistency \
            (current row limit: %d)"
            i rows_i limit ()
        else loop ~i:(i - 1) ~limit:rows_i
    in
    loop ~i:(Slap.Size.to_int n_rows) (*! S2I *) ~limit:real_m
  end

(* Checks whether a sparse column matrix is sane *)
let check_sparse_col_mat_sane ~real_n ~smat ~cols =
  if !debug then begin
    if real_n < 0 then
      failwith "Gpr.Gpr_utils.check_sparse_col_mat_sane: real_n < 0";
    let (n : 'n Slap.Size.t) (*! ITA[4] *) = Mat.dim2 smat in
    let (n_cols : 'n Slap.Size.t) (*! ITA[4] *) = Int_vec.dim cols in
    (* if n_cols <> n then *) (*! RMDC *)
    (*   failwithf *) (*! RMDC *)
    (*     "Gpr.Gpr_utils.check_sparse_col_mat_sane: number of cols in \ *) (*! RMDC *)
    (*     sparse matrix (%d) disagrees with size of col array (%d)" *) (*! RMDC *)
    (*     n n_cols (); *) (*! RMDC *)
    let rec loop ~i ~limit =
      if i > 0 then
        let cols_i = Slap.Vec.get_dyn cols i (*! IDX *) in
        if cols_i <= 0 then
          failwithf
            "Gpr.Gpr_utils.check_sparse_col_mat_sane: sparse col %d contains \
            illegal negative real col index %d" i cols_i ()
        else if cols_i > limit then
          failwithf
            "Gpr.Gpr_utils.check_sparse_col_mat_sane: sparse col %d \
            associated with real col index %d violates consistency \
            (current col limit: %d)"
            i cols_i limit ()
        else loop ~i:(i - 1) ~limit:cols_i
    in
    loop ~i:(Slap.Size.to_int n_cols) (*! S2I *) ~limit:real_n
  end

(* Checks whether a parse vector is sane *)
let check_sparse_vec_sane ~real_n ~(svec : ('n, 'x_cd) vec) ~(rows : ('n, 'y_cd) Int_vec.t) (*! ITA[5] *) =
  (* if !debug then *) (*! RMDC *)
  (*   let k = Vec.dim svec in *) (*! RMDC *)
  (*   if Int_vec.dim rows <> k then *) (*! RMDC *)
  (*     failwith *) (*! RMDC *)
  (*       "Gpr.Gpr_utils.check_sparse_vec_sane: \ *) (*! RMDC *)
  (*       size of sparse vector disagrees with indexes"; *) (*! RMDC *)
    let rec loop ~last i =
      if i > 0 then
        let ind = Slap.Vec.get_dyn rows i (*! IDX *) in
        if ind >= last || ind <= 0 then
          failwith "Gpr.Gpr_utils.check_sparse_vec_sane: rows inconsistent"
        else loop ~last:ind (i - 1)
    in
    loop ~last:real_n (Slap.Size.to_int (Int_vec.dim rows)) (*! S2I *)

(* Computes the trace of the product of a symmetric and sparse
   symmetric matrix *)
let symm2_sparse_trace ~mat ~smat ~rows =
  let m = Slap.Size.to_int (Int_vec.dim rows) (*! S2I *) in
  let n = Slap.Size.to_int (Mat.dim2 smat) (*! S2I *) in
  let full_ref = ref 0. in
  let half_ref = ref 0. in
  let rows_ix_ref = ref 1 in
  for sparse_r = 1 to m do
    let c = Slap.Vec.get_dyn rows sparse_r (*! IDX *) in
    for r = 1 to n do
      let mat_el = if r > c then (Mat.get_dyn mat c r) else (Mat.get_dyn mat r c) in (*! IDX *)
      let rows_ix = !rows_ix_ref in
      if
        rows_ix > m ||
        let rows_el = Slap.Vec.get_dyn rows rows_ix (*! IDX *) in
        r < rows_el || c < rows_el
      then full_ref := !full_ref +. mat_el *. (Mat.get_dyn smat sparse_r r) (*! IDX *)
      else begin
        half_ref := !half_ref +. mat_el *. (Mat.get_dyn smat rows_ix c); (*! IDX *)
        incr rows_ix_ref
      end
    done;
    rows_ix_ref := 1
  done;
  let full = !full_ref in
  full +. !half_ref +. full

type ('D, 'd, 'k, 'cd) id = (('D, 'cd) vec -> ('d, 'cd) vec) * (('D, 'k, 'cd) mat -> ('d, 'k, 'cd) mat) (*! FS[1] *)
