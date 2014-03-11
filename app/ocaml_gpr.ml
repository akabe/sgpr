(* File: ocaml_gpr.ml

   Sized GPR - OCaml-GPR with static size checking of operations on matrices

   [Authors of Sized GPR]
     Copyright (C) 2014-  Akinori ABE
     email: abe@kb.ecei.tohoku.ac.jp

   [Authors of OCaml-GPR]
     Copyright (C) 2009-  Markus Mottl
     email: markus.mottl@gmail.com
     WWW:   http://www.ocaml.info

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*)

open Core.Std

module Args = struct
  type cmd = [ `Train | `Test ]

  type t = {
    cmd : cmd;
    model_file : string;
    with_stddev : bool;
    predictive : bool;
    max_iter : int option;
    n_inducing : int;
    sigma2 : float;
    amplitude : float;
    dim_red : int option;
    log_het_sked : float option;
    multiscale : bool;
    tol : float;
    step : float;
    eps : float;
    verbose : bool;
  }

  let cmd : cmd ref = ref `Train
  let model_file = ref None
  let with_stddev = ref false
  let predictive = ref false
  let max_iter = ref None
  let n_inducing = ref 10
  let sigma2 = ref 1.
  let amplitude = ref 1.
  let dim_red = ref None
  let log_het_sked = ref None
  let multiscale = ref false
  let tol = ref 0.1
  let step = ref 0.1
  let eps = ref 0.1
  let verbose = ref false

  let set_some n_ref n = n_ref := Some n

  let args =
    Arg.align
      [
        (
          "-cmd",
          Arg.Symbol ([ "train"; "test" ], function
            | "train" -> cmd := `Train
            | "test" -> cmd := `Test
            | _ -> assert false  (* impossible *)
          ),
          " train (default) or test model"
        );(
          "-model",
          Arg.String (fun str -> model_file := Some str),
          " model file to use"
        );(
          "-with-stddev",
          Arg.Set with_stddev,
          " make predictions with both mean and variance"
        );(
          "-predictive",
          Arg.Set predictive,
          " standard deviation includes noise level (predictive distribution)"
        );(
          "-max-iter",
          Arg.Int (set_some max_iter),
          " maximum number of optimization steps (default: limitless)"
        );(
          "-n-inducing",
          Arg.Set_int n_inducing,
          sprintf
            " sets number of randomly initialized inducing inputs (default: %d)"
            !n_inducing
        );(
          "-sigma2",
          Arg.Set_float sigma2,
          sprintf " sets initial noise level (default: %f)" !sigma2
        );(
          "-amplitude",
          Arg.Set_float amplitude,
          sprintf " sets initial amplitude level (default: %f)" !amplitude
        );(
          "-dim-red",
          Arg.Int (set_some dim_red),
          " sets dimensionality reduction (default: none)"
        );(
          "-log-het-sked",
          Arg.Float (set_some log_het_sked),
          " turns on / sets log-heteroskedastic \
          noise (may require negative values)"
        );(
          "-multiscale",
          Arg.Set multiscale,
          " turns on multiscale approximation"
        );(
          "-tol",
          Arg.Set_float tol,
          sprintf " sets tolerance for gradient descent (default: %f)" !tol
        );(
          "-step",
          Arg.Set_float step,
          sprintf " sets step size for gradient descent (default: %f)" !step
        );(
          "-eps",
          Arg.Set_float eps,
          sprintf " sets epsilon for gradient descent (default: %f)" !eps
        );(
          "-verbose",
          Arg.Set verbose,
          " prints information while training"
        );
      ]

  let usage_msg = sprintf "%s: -cmd [ train | test ] -model file" Sys.argv.(0)

  let anon_fun _ = failwith "no anonymous arguments allowed"

  let some name opt_ref =
    match !opt_ref with
    | Some v -> v
    | None ->
        eprintf "command line option %s not provided\n\n%!" name;
        prerr_endline usage_msg;
        exit 1

  let get () =
    Arg.parse args anon_fun usage_msg;
    {
      cmd = !cmd;
      model_file = some "model" model_file;
      with_stddev = !with_stddev;
      predictive = !predictive;
      max_iter = !max_iter;
      n_inducing = !n_inducing;
      sigma2 = !sigma2;
      amplitude = !amplitude;
      dim_red = !dim_red;
      log_het_sked = !log_het_sked;
      multiscale = !multiscale;
      tol = !tol;
      step = !step;
      eps = !eps;
      verbose = !verbose;
    }
end

let read_samples () =
  let rex = Str.regexp "," in
  let split str = Array.of_list (Str.split rex str) in
  match try Some (read_line ()) with _ -> None with
  | None -> failwith "no data"
  | Some line ->
      let conv_line line =
        try Array.map ~f:Float.of_string (split line)
        with exc -> Exn.reraisef exc "failure '%s' converting sample" line ()
      in
      let sample = conv_line line in
      let d = Array.length sample in
      let rec loop samples =
        match try Some (read_line ()) with _ -> None with
        | Some line ->
            let floats = conv_line line in
            if Array.length floats <> d then
              failwithf
                "incompatible dimension of sample in line %d: %s"
                (List.length samples + 1) line ()
            else loop (floats :: samples)
        | None -> Array.of_list (List.rev samples)
      in
      loop [sample]

open Slap.D (*! RID *)

open Sgpr (*! RID *)

module GP = Fitc_gp.Make_deriv (Cov_se_fat.Deriv)
module FIC = GP.Variational_FIC.Eval

module Model = struct
  type ('a, 'b, 'c, 'd, 'e, 'f, 'g, 'h, 'i) t = { (*! ITP *)
    sigma2 : float;
    target_mean : float;
    input_means : ('a, Slap.cnt) vec; (*! ITP *)
    input_stddevs : ('b, Slap.cnt) vec; (*! ITP *)
    kernel : ('c, 'd, 'e) Cov_se_fat.Eval.Kernel.t; (*! ITP *)
    inducing_points : ('f, 'g) FIC.Spec.Inducing.t; (*! ITP *)
    coeffs : ('h, Slap.cnt) vec; (*! ITP *)
    co_variance_coeffs : 'i FIC.Model.co_variance_coeffs; (*! ITP *)
  }
end

let read_training_samples samples n d = (*! EGPT[3] *)
(*   let samples = read_samples () in *) (*! EGPT[3] *)
(*   let n = Array.length samples in *) (*! EGPT[3] *)
(*   let d = Array.length samples.(0) - 1 in *) (*! EGPT[3] *)
  let inputs = Mat.create d n in
  let targets = Vec.create n in
  Array.iteri samples ~f:(fun c0 sample ->
    for r1 = 1 to Slap.Size.to_int d (*! S2I *) do Mat.set_dyn inputs r1 (c0 + 1) sample.(r1 - 1) (*! IDX *) done;
    Vec.set_dyn targets (c0 + 1) sample.(Slap.Size.to_int d)); (*! S2I,IDX *)
  inputs, targets

let write_model model_file ~target_mean ~input_means ~input_stddevs trained =
  let oc = open_out model_file in
  let model =
    let model = FIC.Trained.get_model trained in
    let sigma2 = FIC.Model.get_sigma2 model in
    let kernel = FIC.Model.get_kernel model in
    let inducing = FIC.Model.get_inducing model in
    let inducing_points = FIC.Inducing.get_points inducing in
    let mean_predictor = FIC.Mean_predictor.calc_trained trained in
    let coeffs = FIC.Mean_predictor.get_coeffs mean_predictor in
    let co_variance_coeffs = FIC.Model.calc_co_variance_coeffs model in
    {
      Model.
      sigma2; target_mean; input_means; input_stddevs; kernel;
      inducing_points; coeffs; co_variance_coeffs;
    }
  in
  Marshal.to_channel oc model [];
  Out_channel.close oc

exception Bailout

let train args =
  let
    {
      Args.
      model_file; max_iter; n_inducing; sigma2; amplitude; dim_red;
      log_het_sked; multiscale; tol; step; eps = epsabs; verbose
    } = args
  in
  let samples = read_samples () in (*! EGPT[3] *)
  let module N = Slap.Size.Of_int_dyn(struct let value = Array.length samples end) in (*! EGPT[3] *)
  let module D = Slap.Size.Of_int_dyn(struct let value = Array.length samples.(0) - 1 end) in (*! EGPT[3] *)
  let inputs, targets = read_training_samples samples N.value D.value (*! EGPT[3] *) in
  let big_dim = Mat.dim1 inputs in
  let n_inputs = Mat.dim2 inputs in
  let f_inputs = float (Slap.Size.to_int n_inputs) (*! S2I *) in
  let calc_mean vec = Vec.sum vec /. float (Slap.Size.to_int (Vec.dim vec)) (*! S2I *) in
  let target_mean = calc_mean targets in
  let targets = Vec.map (fun n -> n -. target_mean) targets in
  let target_variance = Vec.sqr_nrm2 targets /. f_inputs in
  if verbose then eprintf "target variance: %.5f\n%!" target_variance;
  let input_means = Vec.create big_dim in
  let input_stddevs = Vec.create big_dim in
  for i = 1 to Slap.Size.to_int big_dim (*! S2I *) do
    let input = Mat.copy_row_dyn (*! RID *) inputs i in
    let mean = calc_mean input in
    Vec.set_dyn input_means i mean; (*! IDX *)
    let stddev = sqrt (Vec.ssqr ~c:mean input) in
    Vec.set_dyn input_stddevs i stddev; (*! IDX *)
    for j = 1 to Slap.Size.to_int n_inputs (*! S2I *) do
      Mat.set_dyn inputs i j ((Mat.get_dyn inputs i j -. mean) /. stddev); (*! IDX *)
    done;
  done;
  let module M = Slap.Size.Of_int_dyn(struct let value = n_inducing end) in (*! I2S *)
  let n_inducing = Slap.Size.min M.value (*! I2S *) (Vec.dim targets) (*! SOP *) in
  Random.self_init ();
  let f d tproj = (*! ET[1] *)
    let log_sf2 = 2. *. log amplitude in
    let params = (*! ET[1] *)
      let log_hetero_skedasticity =
        match log_het_sked with
        | Some log_het_sked -> Some (Vec.make n_inducing log_het_sked)
        | None -> None
      in
      let log_multiscales_m05 =
        if multiscale then Some (Mat.make0 d n_inducing)
        else None
      in
      Cov_se_fat.Params.create
        {
          Cov_se_fat.Params.
          d; log_sf2; tproj; log_hetero_skedasticity; log_multiscales_m05
        }
    in
    let kernel = Cov_se_fat.Eval.Kernel.create params in
    let get_trained_stats trained =
      let { FIC.Stats.smse; msll; mad; maxad } = FIC.Stats.calc trained in
      sprintf
        "MSLL=%7.7f SMSE=%7.7f MAD=%7.7f MAXAD=%7.7f"
        msll smse mad maxad
    in
    let best_trained = ref None in
    let report_trained_model, report_gradient_norm =
      let got_signal = ref false in
      Signal.Expert.set Signal.int (`Handle (fun _ -> got_signal := true));
      let bailout ~iter _ =
        if !got_signal then raise Bailout;
        match max_iter with
        | Some max_iter when iter > max_iter -> raise Bailout
        | _ -> ()
      in
      if verbose then
        let last_eval_time = ref 0. in
        let last_deriv_time = ref 0. in
        let maybe_print last_time line =
          let now = Unix.gettimeofday () in
          if !last_time +. 1. < now then begin
                                        last_time := now;
                                        prerr_endline line;
                                      end
        in
        Some (fun ~iter trained ->
              best_trained := Some trained;
              bailout ~iter ();
              maybe_print last_eval_time
                          (sprintf "iter %4d: %s" iter (get_trained_stats trained))),
        Some (fun ~iter norm ->
              bailout ~iter ();
              maybe_print last_deriv_time
                          (sprintf "iter %4d: |gradient|=%.5f" iter norm))
      else Some bailout, None
    in
    let hypers = GP.FIC.Deriv.Optim.get_default_hypers ~kernel ~n_rand_inducing:n_inducing ~inputs in (*! EGPT[2] *)
    let module N = Slap.Size.Of_int_dyn(struct let value = Array.length hypers end) in (*! EGPT[2] *)
    match
      try
        Some (
            GP.Variational_FIC.Deriv.Optim.Gsl.train
              ?report_trained_model ?report_gradient_norm
              ~kernel ~sigma2 ~n_rand_inducing:n_inducing
              ~hypers ~n_hypers:N.value (*! EGPT[2] *) ~tol ~step ~epsabs ~inputs ~targets ())
      with GP.FIC.Deriv.Optim.Gsl.Optim_exception Bailout -> !best_trained
    with
    | None -> ()
    | Some trained ->
       if verbose then eprintf "result: %s\n%!" (get_trained_stats trained);
       write_model model_file ~target_mean ~input_means ~input_stddevs trained
  in (*! ET[1] *)
  match dim_red with (*! ET[1] *)
  | None -> f big_dim Cov_se_fat.Params.tproj_none (*! ET[1] *)
  | Some small_dim -> (*! ET[1] *)
     let module K = Slap.Size.Of_int_dyn(struct let value = small_dim end) in (*! I2S *) (*! ET[1] *)
     let small_dim = Slap.Size.min big_dim K.value (*! I2S,SOP *) in (*! ET[1] *)
     let tproj = Mat.random big_dim small_dim in (*! ET[1] *)
     Mat.scal (1. /. float (Slap.Size.to_int big_dim (*! S2I *))) tproj; (*! ET[1] *)
     f small_dim (Cov_se_fat.Params.tproj_some tproj) (*! ET[1] *)

let read_test_samples big_dim samples n (*! EGPT[3] *) =
  (* let samples = read_samples () in *) (*! EGPT[3] *)
  (* let n = Array.length samples in *) (*! EGPT[3] *)
  (* if n = 0 then Mat.empty *) (*! FT[3] *)
  (* else *) begin (*! FT[3] *)
    let input_dim = Array.length samples.(0) in
    if input_dim <> Slap.Size.to_int big_dim (*! S2I *) then
      failwithf
        "incompatible dimension of inputs (%d), expected %d"
        input_dim (Slap.Size.to_int big_dim (*! S2I *)) ();
    let inputs = Mat.create big_dim n in
    Array.iteri samples ~f:(fun c0 sample ->
      for r1 = 1 to Slap.Size.to_int big_dim (*! S2I *) do Mat.set_dyn inputs r1 (c0 + 1) sample.(r1 - 1) (*! IDX *) done);
    inputs
  end

let read_model model_file : (_,_,_,_,_,_,_,_,_) Model.t (*! ITP *) =
  let ic = open_in model_file in
  let model = Marshal.from_channel ic in
  In_channel.close ic;
  model

let test args =
  let { Args.model_file; with_stddev; predictive } = args in
  let
    {
      Model.
      sigma2; target_mean; input_means; input_stddevs; kernel;
      inducing_points; coeffs; co_variance_coeffs
    } = read_model model_file
  in
  let big_dim = Vec.dim input_stddevs in
  (* let inputs = read_test_samples big_dim in *) (*! ET[2],FT[3] *)
  let f inputs = (*! ET[2] *)
    let n_inputs = Mat.dim2 inputs in
    for i = 1 to Slap.Size.to_int big_dim (*! S2I *) do
      let mean = Vec.get_dyn input_means i in (*! IDX *)
      let stddev = Vec.get_dyn input_stddevs i in (*! IDX *)
      for j = 1 to Slap.Size.to_int n_inputs (*! S2I *) do
        Mat.set_dyn inputs i j ((Mat.get_dyn inputs i j -. mean) /. stddev); (*! IDX *)
      done;
    done;
    let mean_predictor = FIC.Mean_predictor.calc inducing_points ~coeffs in
    let inducing = FIC.Inducing.calc kernel inducing_points in
    let inputs = FIC.Inputs.calc inputs inducing in
    let means = FIC.Means.get (FIC.Means.calc mean_predictor inputs) in
    let renorm_mean mean = mean +. target_mean in
    if with_stddev then
      let co_variance_predictor =
        FIC.Co_variance_predictor.calc kernel inducing_points co_variance_coeffs
      in
      let vars = FIC.Variances.calc co_variance_predictor ~sigma2 inputs in
      let vars = FIC.Variances.get ~predictive vars in
      Vec.iteri (fun i pre_mean ->
                 let mean = renorm_mean pre_mean in
                 printf "%f,%f\n" mean (sqrt (Vec.get_dyn vars i) (*! IDX *))) means
    else Vec.iter (fun mean -> printf "%f\n" (renorm_mean mean)) means
  in (*! ET[2] *)
  let samples = read_samples () in (*! EGPT[3] *)
  let module N = Slap.Size.Of_int_dyn(struct let value = Array.length samples end) in (*! EGPT[3] *)
  if Slap.Size.to_int N.value (*! S2I *) = 0 (*! FT[3] *)
  then f Mat.empty (*! FT[3],ET[2] *)
  else f (read_test_samples big_dim samples N.value (*! EGPT[3] *)) (*! FT[3],ET[2] *)

let main () =
  let args = Args.get () in
  match args.Args.cmd with
  | `Train -> train args
  | `Test -> test args

let () = main ()
