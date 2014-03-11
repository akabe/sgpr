(* File: interfaces.ml

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

open Gpr_utils


(** {6 Representations of (sparse) derivative matrices} *)

(** Representation of indices into sparse matrices *)
module Sparse_indices = Int_vec

(** Derivative representations for both symmetric and unsymmetric matrices.

  - Dense: matrix is dense.
  - Sparse_rows: matrix is zero everywhere except for rows whose
    index is stored in the sparse index argument.  The rows in the
    matrix correspond to the given indices.
  - Const: matrix is constant everywhere.
  - Factor: matrix is the non-derived matrix times the given factor
    (useful with exponential functions).
*)
type ('dm, 'sm, 'n, 'cnt_or_dsc) common_mat_deriv = [ (*! ITP *)
  | `Dense of ('dm, 'n, 'cnt_or_dsc) mat (*! ITP *)
  | `Sparse_rows of ('sm, 'n, 'cnt_or_dsc) mat * ('sm, Slap.cnt) Sparse_indices.t (*! ITP *)
  | `Const of float
  | `Factor of float
]

(** Only general matrices support sparse column representations.

  - Sparse_cols: matrix is zero everywhere except for columns whose
    index is stored in the sparse index argument.  The columns in
    the matrix correspond to the given indices.
*)
type ('dm, 'sm, 'dn, 'sn, 'cnt_or_dsc) mat_deriv = [ (*! ITP *)
  | ('dm, 'sm, 'dn, 'cnt_or_dsc) common_mat_deriv (*! ITP *)
  | `Sparse_cols of ('dm, 'sn, 'cnt_or_dsc) mat * ('sn, Slap.cnt) Sparse_indices.t (*! ITP *)
]

(** Only symmetric (square) matrices support diagonal vectors and
    diagonal constants as derivatives.

  - Diag_vec: matrix is zero everywhere except for the diagonal
    whose values are given in the argument.
  - Diag_const: matrix is zero everywhere except for the diagonal
    whose values are set to the given constant.

  Note that sparse rows do not need to compute or store all elements
  for symmetric matrices.  Entries that have already appeared in
  previous rows by symmetry can be left uninitialized.
*)
type ('sm, 'n, 'cnt_or_dsc) symm_mat_deriv = [ (*! ITP *)
  | ('n, 'sm, 'n, 'cnt_or_dsc) common_mat_deriv (*! ITP *)
  | `Diag_vec of ('n, 'cnt_or_dsc) vec (*! ITP *)
  | `Diag_const of float
]

(** Derivatives of diagonal matrices.

  - Vec: the derivatives of the diagonal given in a dense vector.
  - Sparse_vec: matrix is zero everywhere except at those indices
    along the diagonal that are mentioned in the sparse indices
    argument.  The element associated with such an index is stored
    in the vector argument.
  - Const: the derivative of the diagonal matrix is a constant.
  - Factor: the derivative of the diagonal is the the non-derived
    diagonal matrix times the given factor (useful with exponential
    functions).
*)
type ('n, 'cnt_or_dsc) diag_deriv = [ (*! ITP *)
  | `Vec of ('n, 'cnt_or_dsc) vec (*! ITP *)
  | `Sparse_vec of ('n, 'cnt_or_dsc) vec * ('n, Slap.cnt) Sparse_indices.t (*! ITP *)
  | `Const of float
  | `Factor of float
]

(** Specifications of covariance functions (= kernels) and their derivatives *)
module Specs = struct

  (** Signature of kernels and their parameters *)
  module type Kernel = sig
    type ('D, 'd, 'm) t (*! ITP *)  (** Type of kernel *)
    type ('D, 'd, 'm) params (*! ITP *)  (** Type of kernel parameters *)

    (** [create params] @return kernel given parameters [params]. *)
    val create : ('D, 'd, 'm) params -> ('D, 'd, 'm) t (*! ITP *)

    (** [get_params kernel] @return parameters used to parameterize [kernel]. *)
    val get_params : ('D, 'd, 'm) t -> ('D, 'd, 'm) params (*! ITP *)
  end

  (** Evaluation of covariance functions *)
  module type Eval = sig

    (** Kernel used for evaluation *)
    module Kernel : Kernel

    (** Signature for evaluating inducing inputs *)
    module Inducing : sig
      type ('m, 'n) t (*! ITP *)

      (** [get_n_points inducing] @return number of inducing points. *)
      val get_n_points : ('m, 'n) t -> 'n Slap.Size.t (*! ITP *)

      (** [calc_upper kernel inducing] @return upper triangle of
          covariance matrix of [inducing] inputs given [kernel]. *)
      val calc_upper : (_, _, _) Kernel.t -> ('m, 'n) t -> ('n, 'n, 'cnt) mat (*! ITP *)
    end

    (** Signature for evaluating single inputs *)
    module Input : sig
      type 'n t (*! ITP *)  (** Type of input point *)

      (** [eval kernel input inducing] @return (row) vector of covariance
          evaluations between [input] and [inducing] inputs given
          [kernel]. *)
      val eval : ('D, 'd, _) Kernel.t -> 'D t -> ('d, 'n) Inducing.t -> ('n, 'cnt) vec (*! ITP *)

      (** [weighted_eval kernel input inducing ~coeffs] @return
          [coeff]-weighted sum of covariances between [input] and
          [inducing] inputs given [kernel]. *)
      val weighted_eval : ('D, 'd, _) Kernel.t -> 'D t -> ('d, 'n) Inducing.t -> coeffs : ('n, 'cd) vec -> float (*! ITP *)

      (** [eval_one kernel point] @return variance of [point] given [kernel]. *)
      val eval_one : (_, _, _) Kernel.t -> _ t -> float (*! ITP *)
    end

    (** Signature for evaluating multiple inputs *)
    module Inputs : sig
      type ('m, 'n) t (*! ITP *)  (** Type of input points *)

      (** [create inputs] @return inputs given an array of single [inputs]. *)
      val create : 'm Slap.Size.t -> 'n Slap.Size.t (*! EGPT[1] *) -> 'm Input.t array -> ('m, 'n) t (*! ITP *)

      (** [get_n_points inputs] @return number of input points. *)
      val get_n_points : ('m, 'n) t -> 'n Slap.Size.t (*! ITP *)

      (** [choose_subset inputs indexes] @return subset of input
          points from [inputs] having [indexes]. *)
      val choose_subset : ('m, 'n) t -> ('k, 'cd) Int_vec.t -> ('m, 'k) t (*! ITP *)

      (** [create_inducing kernel inputs] @return inducing points
          made from [inputs] and given [kernel]. *)
      val create_inducing : ('D, 'd, _) Kernel.t -> ('D, 'n) t -> ('d, 'n) Inducing.t (*! ITP *)

      type 'D default_kernel_size (*! DKS[1] *)

      (** [create_default_kernel_params inputs ~n_inducing] @return
          default kernel parameters to be used with [n_inducing]
          inducing points and [inputs]. *)
      val create_default_kernel_params : ('D, _) t -> n_inducing : 'm Slap.Size.t -> ('D, 'D default_kernel_size (*! DKS[1] *), 'm) Kernel.params (*! ITP *)

      (** [calc_upper kernel inputs] @return upper triangle of
          covariance matrix of [inputs] given [kernel]. *)
      val calc_upper : ('D, _, _) Kernel.t -> ('D, 'n) t -> ('n, 'n, 'cnt) mat (*! ITP *)

      (** [calc_diag kernel inputs] @return diagonal of
          covariance matrix of [inputs] given [kernel]. *)
      val calc_diag : ('D, _, _) Kernel.t -> ('D, 'n) t -> ('n, 'cnt) vec (*! ITP *)

      (** [calc_cross kernel ~inputs ~inducing] @return cross-covariance
          matrix of [inputs] (indexing rows) and [inducing] points
          (indexing columns). *)
      val calc_cross : ('D, 'd, _) Kernel.t -> inputs : ('D, 'm) t -> inducing : ('d, 'n) Inducing.t -> ('m, 'n, 'cnt) mat (*! ITP *)

      (** [weighted_eval kernel ~inputs ~inducing ~coeffs] @return
          vector of [coeff]-weighted sums of covariances between
          [inputs] and [inducing] inputs given [kernel]. *)
      val weighted_eval :
        ('D, 'd, _) Kernel.t -> inputs : ('D, 'n) t -> inducing : ('d, 'n) Inducing.t -> coeffs : ('n, 'cd) vec -> ('n, 'cnt) vec (*! ITP *)
    end
  end

  (** Derivatives of covariance functions *)
  module type Deriv = sig

    (** Derivatives always require evaluation functions *)
    module Eval : Eval

    (** Hyper parameters that have derivatives *)
    module Hyper : sig
      type t  (** Type of hyper parameter *)

      (** [get_all kernel inducing inputs] @return array of all hyper
          parameters of [kernel] and/or ([inducing]) [inputs] for which
          derivatives can be computed. *)
      val get_all : (_, _, _) Eval.Kernel.t -> (_, _) Eval.Inducing.t -> (_, _) Eval.Inputs.t -> t array (*! ITP *)

      (** [get_value kernel inducing inputs hyper] @return value of hyper
          parameter [hyper] of [kernel] and/or ([inducing]) [inputs]. *)
      val get_value :
        (_, _, _) Eval.Kernel.t -> (_, _) Eval.Inducing.t -> (_, _) Eval.Inputs.t -> t -> float (*! ITP *)

      (** [set_values kernel inducing inputs hypers values] @return triple
          of [(kernel, inducing, inputs)] in which [hypers] have been
          substituted with [values] position-wise. *)
      val set_values :
        ('D, 'd, 'm) Eval.Kernel.t -> ('idm, 'idn) Eval.Inducing.t -> ('ipm, 'ipn) Eval.Inputs.t -> t array -> (_, _) vec -> (*! ITP *)
        ('D, 'd, 'm) Eval.Kernel.t * ('idm, 'idn) Eval.Inducing.t * ('ipm, 'ipn) Eval.Inputs.t (*! ITP *)
    end

    (** Derivatives of the covariance matrix of inducing inputs *)
    module Inducing : sig
      type ('D, 'd, 'm, 'n) upper (*! ITP *)  (** Representation of precomputed data for calculating the
                      upper triangle of the derivative of the covariance matrix
                      of inducing inputs. *)

      (** [calc_shared_upper kernel inducing] @return the pair [(eval, upper)],
          where [eval] is the upper triangle of the covariance matrix of
          inducing inputs for [kernel], and [upper] is the precomputed data
          needed for taking derivatives. *)
      val calc_shared_upper : ('D, 'd, 'm) Eval.Kernel.t -> ('n, 'm) Eval.Inducing.t -> ('m, 'm, Slap.cnt) mat * ('D, 'd, 'm, 'n) upper (*! ITP *)

      (** [calc_deriv_upper upper hyper] @return the derivative of the
          (symmetric) covariance matrix of inducing inputs given precomputed
          data [upper] and the [hyper]-variable. *)
      val calc_deriv_upper : ('D, 'd, 'm, 'n) upper -> Hyper.t -> (Slap.Size.z Slap.Size.s, 'm, 'cnt) symm_mat_deriv (*! ITP *)
    end

    (** Derivatives of the (cross-) covariance matrix of inputs. *)
    module Inputs : sig
      type ('D, 'd, 'm, 'n) diag (*! ITP *)  (** Representation of precomputed data for calculating the
                     derivative of the diagonal of the covariance matrix of
                     inputs. *)

      type ('D, 'd, 'm, 'n) cross (*! ITP *)  (** Representation of precomputed data for calculating the
                      derivative of the cross-covariance matrix between inputs
                      and inducing inputs. *)

      (** [calc_shared_diag kernel inputs] @return the pair [(eval, diag)],
          where [eval] is the diagonal of the covariance matrix of [inputs] for
          [kernel], and [diag] is the precomputed data needed for taking
          derivatives. *)
      val calc_shared_diag : ('D, 'd, 'm) Eval.Kernel.t -> ('D, 'n) Eval.Inputs.t -> ('n, 'cnt) vec * ('D, 'd, 'm, 'n) diag (*! ITP *)

      (** [calc_shared_cross kernel ~inputs ~inducing] @return the pair [(eval,
          cross)], where [eval] is the cross-covariance matrix of inputs and
          inducing inputs for [kernel], and [diag] is the precomputed data
          needed for taking derivatives. *)
      val calc_shared_cross :
        ('D, 'd, 'm) Eval.Kernel.t -> inputs : ('D, 'n) Eval.Inputs.t -> inducing : ('d, 'm) Eval.Inducing.t -> (*! ITP *)
        ('n, 'm, Slap.cnt) mat * ('D, 'd, 'm, 'n) cross (*! ITP *)

      (** [calc_deriv_diag diag hyper] @return the derivative of the
          diagonal of the covariance matrix of inputs given precomputed data
          [diag] and the [hyper]-variable. *)
      val calc_deriv_diag : (_, _, _, 'n) diag -> Hyper.t -> ('n, 'cnt) diag_deriv (*! ITP *)

      (** [calc_deriv_cross cross hyper] @return the derivative of the
          cross-covariance matrix of the inputs and inducing inputs given
          precomputed data [cross] and the [hyper]-variable. *)
      val calc_deriv_cross : ('D, 'd, 'm, 'n) cross -> Hyper.t -> ('n, 'sm, 'm, Slap.Size.z Slap.Size.s, 'cnt) mat_deriv (*! ITP *)
    end
  end

  (** Derivatives of inputs for global optimization. *)
  module type Optimizer = sig

    (** Derivatives always require evaluation functions *)
    module Eval : Eval

    (** Input parameters that have derivatives *)
    module Var : sig
      type t  (** Type of input parameter *)
    end

    module Input : sig
      (** [get_vars input] @return array of all input parameters for which
          derivatives can be computed given [input]. *)
      val get_vars : _ Eval.Input.t -> Var.t array (*! ITP *)

      (** [get_value input var] @return value of input parameter [var] for
          [input]. *)
      val get_value : _ Eval.Input.t -> Var.t -> float (*! ITP *)

      (** [set_values input vars values] @return input in which [vars] have been
          substituted with [values] position-wise. *)
      val set_values : _ Eval.Input.t -> Var.t array -> (_, _) vec -> _ Eval.Input.t (*! ITP *)
    end

    module Inputs : sig
      (** [get_vars inputs] @return array of all input parameters for which
          derivatives can be computed given [inputs]. *)
      val get_vars : (_, _) Eval.Inputs.t -> Var.t array (*! ITP *)

      (** [get_value inputs var] @return value of input parameter [var] for
          [inputs]. *)
      val get_value : (_, _) Eval.Inputs.t -> Var.t -> float (*! ITP *)

      (** [set_values inputs vars values] @return inputs in which [vars] have
          been substituted with [values] position-wise. *)
      val set_values : (_, _) Eval.Inputs.t -> Var.t array -> (_, _) vec -> (_, _) Eval.Inputs.t (*! ITP *)
    end
  end
end

(** Signatures for learning sparse Gaussian processes with inducing inputs *)
module Sigs = struct

  (** Modules for learning without derivatives of covariance functions. *)
  module type Eval = sig

    (** Specification of covariance function *)
    module Spec : Specs.Eval

    (** Evaluating inducing inputs *)
    module Inducing : sig
      type ('D, 'd, 'm) t (*! ITP *)  (** Type of inducing inputs *)

      (** [choose_n_first_inputs kernel inputs ~n_inducing] @return the first
          [n_inducing] inputs in [inputs] as inducing points given [kernel]. *)
      val choose_n_first_inputs :
        ('D, 'd, _) Spec.Kernel.t -> ('D, _) Spec.Inputs.t -> n_inducing : 'm Slap.Size.t -> ('d, 'm) Spec.Inducing.t (*! ITP *)

      (** [choose_n_random_inputs ?rnd_state kernel inputs ~n_inducing] @return
          [n_inducing] random inputs in [inputs] as inducing points given
          [kernel] and (optional) random state [rnd_state].

          @param rnd_state default = default used by the Random module
      *)
      val choose_n_random_inputs :
        ?rnd_state : Random.State.t ->
        ('D, 'd, _) Spec.Kernel.t -> (*! ITP *)
        ('D, _) Spec.Inputs.t -> (*! ITP *)
        n_inducing : 'n Slap.Size.t -> (*! ITP *)
        ('d, 'n) Spec.Inducing.t (*! ITP *)

      (** [calc kernel inducing_points] @return inducing inputs (= precomputed
          data) prepared using [inducing_points] and [kernel]. *)
      val calc : ('D, 'd, 'm) Spec.Kernel.t -> ('d, 'm) Spec.Inducing.t -> ('D, 'd, 'm) t (*! ITP *)

      (** [get_points kernel inducing] @return inducing points associated with
          the prepared [inducing] inputs. *)
      val get_points : ('D, 'd, 'm) t -> ('d, 'm) Spec.Inducing.t (*! ITP *)
    end

    (** Evaluating single inputs *)
    module Input : sig
      type ('D, 'd, 'm) t (*! ITP *)  (** Type of single input *)

      (** [calc inducing point] @return input (= precomputed
          data) prepared using [inducing] inputs and input [point]. *)
      val calc : ('D, 'd, 'm) Inducing.t -> 'D Spec.Input.t -> ('D, 'd, 'm) t (*! ITP *)
    end

    (** Evaluating (multiple) inputs *)
    module Inputs : sig
      type ('D, 'd, 'm, 'n) t (*! ITP *)  (** Type of (multiple) inputs *)

      (** [create_default_kernel points] @return a default kernel given input
          [points] and [n_inducing] inducing inputs. *)
      val create_default_kernel :
        ('D, _) Spec.Inputs.t -> n_inducing : 'n Slap.Size.t -> ('D, 'D Spec.Inputs.default_kernel_size (*! DKS[1] *), 'n) Spec.Kernel.t (*! ITP *)

      (** [create points inducing] @return inputs (= precomputed
          data) prepared using [inducing] inputs and input [points]. *)
      val calc : ('D, 'n) Spec.Inputs.t -> ('D, 'd, 'm) Inducing.t -> ('D, 'd, 'm, 'n) t (*! ITP *)

      (** [get_points kernel inputs] @return points associated with
          the prepared [inputs]. *)
      val get_points : ('D, 'd, 'm, 'n) t -> ('D, 'n) Spec.Inputs.t (*! ITP *)
    end

    (** (Untrained) model - does not require targets *)
    module Model : sig
      type ('D, 'd, 'm, 'n) t (*! ITP *)  (** Type of models *)

      type 'n co_variance_coeffs (*! ITP *)  (** Type of covariance coefficients *)

      (** [calc inputs ~sigma2] @return model given [inputs] and noise level
          [sigma2] (= variance, i.e. squared standard deviation). *)
      val calc : ('D, 'd, 'm, 'n) Inputs.t -> sigma2 : float -> ('D, 'd, 'm, 'n) t (*! ITP *)

      (** [update_sigma2 model sigma2] @return model by updating [model] with
          new noise level [sigma2]. *)
      val update_sigma2 : ('D, 'd, 'm, 'n) t -> float -> ('D, 'd, 'm, 'n) t (*! ITP *)

      (** [calc_log_evidence model] @return the contribution to the log evidence
          (= log marginal likelihood) of [model]. *)
      val calc_log_evidence : ('D, 'd, 'm, 'n) t -> float (*! ITP *)

      (** [calc_co_variance_coeffs model] @return the coefficients
          required for computing posterior (co-)variances for [model]. *)
      val calc_co_variance_coeffs : ('D, 'd, 'm, 'n) t -> 'm co_variance_coeffs (*! ITP *)

      (** [get_kernel model] @return the kernel associated with [model]. *)
      val get_kernel : ('D, 'd, 'm, _) t -> ('D, 'd, 'm) Spec.Kernel.t (*! ITP *)

      (** [get_sigma2 model] @return the noise level associated with [model]. *)
      val get_sigma2 : ('D, 'd, 'm, 'n) t -> float (*! ITP *)

      (** [get_inputs model] @return the inputs associated with [model]. *)
      val get_inputs : ('D, 'd, 'm, 'n) t -> ('D, 'd, 'm, 'n) Inputs.t (*! ITP *)

      (** [get_inputs model] @return the inducing inputs associated with
          [model]. *)
      val get_inducing : ('D, 'd, 'm, 'n) t -> ('D, 'd, 'm) Inducing.t (*! ITP *)
    end

    (** Trained model - requires targets *)
    module Trained : sig
      type ('D, 'd, 'm, 'n) t (*! ITP *)  (** Type of trained models *)

      (** [calc model ~targets] @return trained model given [model] and
          [targets]. *)
      val calc : ('D, 'd, 'm, 'n) Model.t -> targets : ('n, Slap.cnt) vec -> ('D, 'd, 'm, 'n) t (*! ITP *)

      (** [calc_mean_coeffs trained] @return the vector of coefficients for
          computing posterior means. *)
      val calc_mean_coeffs : ('D, 'd, 'm, 'n) t -> ('m, Slap.cnt) vec (*! ITP *)

      (** [calc_log_evidence trained] @return the log evidence for the trained
          model (includes contribution to log evidence by underlying model). *)
      val calc_log_evidence : ('D, 'd, 'm, 'n) t -> float (*! ITP *)

      (** [get_model trained] @return the model associated with the [trained]
          model. *)
      val get_model : ('D, 'd, 'm, 'n) t -> ('D, 'd, 'm, 'n) Model.t (*! ITP *)

      (** [get_targets trained] @return targets used for training [trained]. *)
      val get_targets : ('D, 'd, 'm, 'n) t -> ('n, Slap.cnt) vec (*! ITP *)
    end

    (** Statistics derived from trained models *)
    module Stats : sig
      (** Type of full statistics *)
      type t = {
        n_samples : int;  (** Number of samples used for training *)
        target_variance : float;  (** Variance of targets *)
        sse : float;   (** Sum of squared errors *)
        mse : float;  (** Mean sum of squared errors *)
        rmse : float;  (** Root mean sum of squared errors *)
        smse : float;  (** Standardized mean squared error *)
        msll : float;  (** Mean standardized log loss *)
        mad : float;  (** Mean absolute deviation *)
        maxad : float;  (** Maximum absolute deviation *)
      }

      (** [calc_n_samples trained] @return number of samples used for training
          [trained]. *)
      val calc_n_samples : ('D, 'd, 'm, 'n) Trained.t -> 'n Slap.Size.t (*! ITP *)

      (** [calc_target_variance trained] @return variance of targets used for
          training [trained]. *)
      val calc_target_variance : ('D, 'd, 'm, 'n) Trained.t -> float (*! ITP *)

      (** [calc_sse trained] @return the sum of squared errors of the [trained]
          model. *)
      val calc_sse : ('D, 'd, 'm, 'n) Trained.t -> float (*! ITP *)

      (** [calc_mse trained] @return the mean sum of squared errors of the
          [trained] model. *)
      val calc_mse : ('D, 'd, 'm, 'n) Trained.t -> float (*! ITP *)

      (** [calc_sse trained] @return the root of the mean sum of squared errors
          of the [trained] model. *)
      val calc_rmse : ('D, 'd, 'm, 'n) Trained.t -> float (*! ITP *)

      (** [calc_smse trained] @return the standardized mean squared error of the
          [trained] model.  This is equivalent to the mean squared error divided
          by the target variance. *)
      val calc_smse : ('D, 'd, 'm, 'n) Trained.t -> float (*! ITP *)

      (** [calc_msll trained] @return the mean standardized log loss.  This
          is equivalent to subtracting the log evidence of the trained model
          from the log evidence of a normal distribution fit to the targets, and
          dividing the result by the number of samples. *)
      val calc_msll : ('D, 'd, 'm, 'n) Trained.t -> float (*! ITP *)

      (** [calc_mad trained] @return the mean absolute deviation
          of the [trained] model. *)
      val calc_mad : ('D, 'd, 'm, 'n) Trained.t -> float (*! ITP *)

      (** [calc_mad trained] @return the maximum absolute deviation
          of the [trained] model. *)
      val calc_maxad : ('D, 'd, 'm, 'n) Trained.t -> float (*! ITP *)

      (** [calc trained] @return the full set of statistics associated with
          the [trained] model. *)
      val calc : ('D, 'd, 'm, 'n) Trained.t -> t (*! ITP *)
    end

    (** Module for making mean predictions *)
    module Mean_predictor : sig
      type ('m, 'n) t (*! ITP *)  (** Type of mean predictors *)

      (** [calc inducing_points ~coeffs] @return a mean predictor given
          [inducing_points] and coefficients [coeffs]. *)
      val calc : ('m, 'n) Spec.Inducing.t -> coeffs : ('n, Slap.cnt) vec -> ('m, 'n) t (*! ITP *)

      (** [calc_trained trained] @return a mean predictor given the [trained]
          model. *)
      val calc_trained : ('D, 'd, 'm, 'n) Trained.t -> ('d, 'm) t (*! ITP *)

      (** [get_inducing mean_predictor] @return inducing points associated with
          [mean_predictor]. *)
      val get_inducing : ('m, 'n) t -> ('m, 'n) Spec.Inducing.t (*! ITP *)

      (** [get_coeffs mean_predictor] @return coefficients associated with
          [mean_predictor]. *)
      val get_coeffs : ('m, 'n) t -> ('n, Slap.cnt) vec (*! ITP *)
    end

    (** Posterior mean for a single input *)
    module Mean : sig
      type 'n t (*! ITP *)  (** Type of mean *)

      (** [calc mean_predictor input] @return mean for [input] given
          [mean_predictor]. *)
      val calc : ('d, 'm) Mean_predictor.t -> ('D, 'd, 'm) Input.t -> 'D t (*! ITP *)

      (** [get mean] @return the mean as a float. *)
      val get : 'n t -> float (*! ITP *)
    end

    (** Posterior means for (multiple) inputs *)
    module Means : sig
      type ('m, 'n) t (*! ITP *)  (** Type of means *)

      (** [calc mean_predictor inputs] @return means for [inputs] given
          [mean_predictor]. *)
      val calc : ('d, 'm) Mean_predictor.t -> ('D, 'd, 'm, 'n) Inputs.t -> ('D, 'n) t (*! ITP *)

      (** [get means] @return the means as a vector. *)
      val get : ('m, 'n) t -> ('n, Slap.cnt) vec (*! ITP *)
    end

    (** Module for making (co-)variance predictions *)
    module Co_variance_predictor : sig
      type ('D, 'd, 'm, 'n) t (*! ITP *)  (** Type of (co-)variance predictor *)

      (** [calc kernel inducing_points co_variance_coeffs] @return (co-)variance
          predictor given [kernel], [inducing_points], and the (co-)variance
          coefficients [co_variance_coeffs]. *)
      val calc :
        ('D, 'd, 'm) Spec.Kernel.t -> ('n, 'm) Spec.Inducing.t -> 'm Model.co_variance_coeffs ->  ('D, 'd, 'm, 'n) t (*! ITP *)

      (** [calc_model model] @return (co-)variance predictor given the
          (untrained) [model]. *)
      val calc_model : ('D, 'd, 'm, 'n) Model.t -> ('D, 'd, 'm, 'd) t (*! ITP *)
    end

    (** Posterior variance for a single input *)
    module Variance : sig
      type 'n t (*! ITP *)  (** Type of variance *)

      (** [calc co_variance_predictor ~sigma2 input] @return variance for
          [input] given [mean_predictor] and noise level [sigma2]. *)
      val calc : (_, _, 'm, 'd) Co_variance_predictor.t -> sigma2 : float -> ('D, 'd, 'm) Input.t -> 'D t (*! ITP *)

      (** [get ?predictive variance] @return the [variance] as a float.
          If [predictive] is [true], then the noise level will be added.

          @param predictive default = [true]
      *)
      val get : ?predictive : bool -> 'n t -> float (*! ITP *)
    end

    (** Posterior variances for (multiple) inputs *)
    module Variances : sig
      type ('m, 'n) t (*! ITP *)  (** Type of variances *)

      (** [calc_model_inputs model] @return variances for all inputs used in
          [model]. *)
      val calc_model_inputs : ('D, 'd, 'm, 'n) Model.t -> ('D, 'n) t (*! ITP *)

      (** [calc co_variance_predictor ~sigma2 inputs] @return variances for
          [inputs] given [co_variance_predictor] and noise level [sigma2]. *)
      val calc : (_, _, 'm, 'd) Co_variance_predictor.t -> sigma2 : float -> ('D, 'd, 'm, 'n) Inputs.t -> ('D, 'n) t (*! ITP *)

      (** [get ?predictive variances] @return the [variances] as a vector.
          If [predictive] is [true], then the noise level will be added.

          @param predictive default = [true]
      *)
      val get : ?predictive : bool -> ('m, 'n) t -> ('n, Slap.cnt) vec (*! ITP *)
    end

    (** Posterior covariances *)
    module Covariances : sig
      type ('m, 'n) t (*! ITP *)  (** Type of covariances *)

      (** [calc_model_inputs model] @return covariances for all inputs used in
          [model].  This may be extremely expensive (O(N^2)) for large numbers
          of model inputs. *)
      val calc_model_inputs : ('D, 'd, 'm, 'n) Model.t -> ('D, 'n) t (*! ITP *)

      (** [calc co_variance_predictor ~sigma2 inputs] @return posterior
          covariances for [inputs] given [co_variance_predictor] and noise level
          [sigma2].  This may be extremely expensive (O(N^2)) for large numbers
          of inputs. *)
      val calc : (_, _, 'm, 'd) Co_variance_predictor.t -> sigma2 : float -> ('D, 'd, 'm, 'n) Inputs.t -> ('D, 'n) t (*! ITP *)

      (** [get ?predictive covariances] @return the [covariances] as a matrix.
          If [predictive] is [true], then the noise level will be added (to the
          diagonal only).

          @param predictive default = [true]
      *)
      val get : ?predictive : bool -> ('m, 'n) t -> ('n, 'n, Slap.cnt) mat (*! ITP *)

      (** [get_variances covariances] @return the variances in [covariances]. *)
      val get_variances : ('m, 'n) t -> ('m, 'n) Variances.t (*! ITP *)
    end

    (** Module for sampling single points from the posterior distribution *)
    module Sampler : sig
      type t  (** Type of sampler *)

      (** [calc ?predictive mean variance] @return sampler given [mean] and
          [variance].  If [predictive] is true, the samples will be noisy. *)
      val calc : ?predictive : bool -> 'n Mean.t -> 'n Variance.t -> t (*! ITP *)

      (** [sample ?rng sampler] @return a sample from the posterior distribution
          given [sampler] and GSL random number generator [rng].

          @param rng default = GSL default
      *)
      val sample : ?rng : Gsl.Rng.t -> t -> float

      (** [samples ?rng sampler ~n] @return [n] samples from the posterior
          distribution given [sampler].

          @param rng default = GSL default
      *)
      val samples : ?rng : Gsl.Rng.t -> t -> n : 'n Slap.Size.t -> ('n, 'cnt) vec (*! ITP *)
    end

    (** Module for sampling (multiple) points from the posterior distribution
        accounting for their covariance *)
    module Cov_sampler : sig
      type 'n t (*! ITP *)  (** Type of covariance sampler *)

      (** [calc ?predictive mean variance] @return sampler given [means] and
          [covariances].  If [predictive] is true, the samples will be noisy. *)
      val calc : ?predictive : bool -> ('m, 'n) Means.t -> ('m, 'n) Covariances.t -> 'n t (*! ITP *)

      (** [sample ?rng sampler] @return a sample vector from the posterior
          distribution given [sampler] and GSL random number generator [rng].

          @param rng default = GSL default
      *)
      val sample : ?rng : Gsl.Rng.t -> 'n t -> ('n, 'cnt) vec (*! ITP *)

      (** [samples ?rng sampler ~n] @return matrix of [n] sample vectors (stored
          row-wise) from the posterior distribution given [sampler].

          @param rng default = GSL default
      *)
      val samples : ?rng : Gsl.Rng.t -> 'm t -> n : 'n Slap.Size.t -> ('m, 'n, 'cd) mat (*! ITP *)
    end
  end

  (** Modules for learning with derivatives of the log evidence (evidence
      maximization framework) *)
  module type Deriv = sig

    (** Sub-modules for learning without derivatives. *)
    module Eval : Eval

    (** Sub-modules for learning with derivatives. *)
    module Deriv : sig

      (** Specification of covariance function derivatives *)
      module Spec : Specs.Deriv with module Eval = Eval.Spec

      (** Module for inducing inputs with derivatives *)
      module Inducing : sig
        type ('D, 'd, 'm) t (*! ITP *)  (** Type of inducing inputs with derivatives *)

        (** [calc kernel inducing_points] @return inducing inputs with
            derivative information given [kernel] and [inducing_points]. *)
        val calc : ('D, 'd, 'm) Eval.Spec.Kernel.t -> ('d, 'm) Eval.Spec.Inducing.t -> ('D, 'd, 'm) t (*! ITP *)

        (** [calc_eval inducing] @return inducing inputs without derivative
            information. *)
        val calc_eval : ('D, 'd, 'm) t -> ('D, 'd, 'm) Eval.Inducing.t (*! ITP *)
      end

      (** Module for inputs with derivatives *)
      module Inputs : sig
        type ('D, 'd, 'm, 'n) t (*! ITP *)  (** Type of inputs with derivatives *)

        (** [calc inducing points] @return inputs with derivative information
            given [inducing] inputs and input [points]. *)
        val calc : ('D, 'd, 'm) Inducing.t -> ('D, 'n) Eval.Spec.Inputs.t -> ('D, 'd, 'm, 'n) t (*! ITP *)

        (** [calc_eval inputs] @return inputs without derivative information. *)
        val calc_eval : ('D, 'd, 'm, 'n) t -> ('D, 'd, 'm, 'n) Eval.Inputs.t (*! ITP *)
      end

      (** (Untrained) model with derivative information *)
      module Model : sig
        type ('D, 'd, 'm, 'n) t (*! ITP *)  (** Type of models with derivatives *)
        type ('D, 'd, 'm, 'n) hyper_t (*! ITP *)  (** Type of models for general hyper parameters *)

        (** [calc inputs ~sigma2] @return model with derivative information
            given [inputs] and noise level [sigma2]. *)
        val calc : ('D, 'd, 'm, 'n) Inputs.t -> sigma2 : float -> ('D, 'd, 'm, 'n) t (*! ITP *)

        (** [update_sigma2 model sigma2] @return model with derivative
            information by updating [model] with new noise level [sigma2]. *)
        val update_sigma2 : ('D, 'd, 'm, 'n) t -> float -> ('D, 'd, 'm, 'n) t (*! ITP *)

        (** [calc_eval model] @return model without derivative information given
            [model]. *)
        val calc_eval : ('D, 'd, 'm, 'n) t -> ('D, 'd, 'm, 'n) Eval.Model.t (*! ITP *)

        (** [calc_log_evidence_sigma2 model] @return the derivative of the
            log evidence of [model] with respect to the noise level (sigma2). *)
        val calc_log_evidence_sigma2 : ('D, 'd, 'm, 'n) t -> float (*! ITP *)

        (** [prepare_hyper model] @return the model prepared for calculating
            derivatives for arbitrary hyper parameters. *)
        val prepare_hyper : ('D, 'd, 'm, 'n) t -> ('D, 'd, 'm, 'n) hyper_t (*! ITP *)

        (** [calc_log_evidence hyper_t hyper] @return the derivative of the log
            evidence given prepared model [hyper_t] with respect to the [hyper]
            variable. *)
        val calc_log_evidence : ('D, 'd, 'm, 'n) hyper_t -> Spec.Hyper.t -> float (*! ITP *)
      end

      (** Trained model with derivative information *)
      module Trained : sig
        type ('D, 'd, 'm, 'n) t (*! ITP *)  (** Type of trained models with derivatives *)
        type ('D, 'd, 'm, 'n) hyper_t (*! ITP *)  (** Type of trained models for general hyper parameters *)

        (** [calc model ~targets] @return trained model with derivative
            information given the untrained [model] and [targets]. *)
        val calc : ('D, 'd, 'm, 'n) Model.t -> targets : ('n, Slap.cnt) vec -> ('D, 'd, 'm, 'n) t (*! ITP *)

        (** [calc_eval trained] @return trained model without derivative
            information given [trained]. *)
        val calc_eval : ('D, 'd, 'm, 'n) t -> ('D, 'd, 'm, 'n) Eval.Trained.t (*! ITP *)

        (** [calc_log_evidence_sigma2 trained] @return the derivative of the
            log evidence for the [trained] model with respect to the noise level
            (sigma2).  This includes the contribution to the derivative by
            [model]. *)
        val calc_log_evidence_sigma2 : ('D, 'd, 'm, 'n) t -> float (*! ITP *)

        (** [prepare_hyper trained] @return the trained model prepared for
            calculating derivatives for arbitrary hyper parameters. *)
        val prepare_hyper : ('D, 'd, 'm, 'n) t -> ('D, 'd, 'm, 'n) hyper_t (*! ITP *)

        (** [calc_log_evidence hyper_t hyper] @return the derivative of the log
            evidence given prepared, trained model [hyper_t] with respect to the
            [hyper] variable. *)
        val calc_log_evidence : ('D, 'd, 'm, 'n) hyper_t -> Spec.Hyper.t -> float (*! ITP *)
      end

      (** Module for testing derivative code *)
      module Test : sig
        (** [check_deriv_hyper ?eps ?tol kernel inducing_points points hyper]
            will raise [Failure] if the derivative code provided in the
            specification of the covariance function given parameter [hyper],
            the [kernel], [inducing_points] and input [points] exceeds the
            tolerance [tol] when compared to finite differences using epsilon
            [eps].  The failure exception will contain details on which
            derivative matrix was incorrect and indicate the matrix element.

            @param eps default = [1e-8]
            @param tol default = [1e-2]
        *)
        val check_deriv_hyper :
          ?eps : float ->
          ?tol : float ->
          ('D, 'd, 'm) Eval.Spec.Kernel.t -> (*! ITP *)
          ('d, 'm) Eval.Spec.Inducing.t -> (*! ITP *)
          ('D, 'n) Eval.Spec.Inputs.t -> (*! ITP *)
          Spec.Hyper.t ->
          unit

        (** [self_test ?eps ?tol kernel inducing_points points ~sigma2 ~targets
            hyper] will raise [Failure] if the internal derivative code for the
            log evidence given parameter [hyper], the [kernel],
            [inducing_points], input [points], noise level [sigma2] and
            [targets] exceeds the tolerance [tol] when compared to finite
            differences using epsilon [eps].

            @param eps default = [1e-8]
            @param tol default = [1e-2]
        *)
        val self_test :
          ?eps : float ->
          ?tol : float ->
          ('D, 'd, 'm) Eval.Spec.Kernel.t -> (*! ITP *)
          ('d, 'm) Eval.Spec.Inducing.t -> (*! ITP *)
          ('D, 'n) Eval.Spec.Inputs.t -> (*! ITP *)
          sigma2 : float ->
          targets : ('n, Slap.cnt) vec -> (*! ITP *)
          [ `Sigma2 | `Hyper of Spec.Hyper.t ] ->
          unit
      end

      (** Optimization module for evidence maximization *)
      module Optim : sig

        type 'n default_n_rand_inducing (*! O2L[1] *)
        val default_n_rand_inducing : 'n Slap.Size.t -> 'n default_n_rand_inducing Slap.Size.t (*! O2L[1] *)

        val get_default_hypers : (*! EGPT[2] *)
          kernel:('D, 'd, _) Eval.Spec.Kernel.t -> (*! EGPT[2] *)
          n_rand_inducing:_ Slap.Size.t -> (*! EGPT[2] *)
          inputs:('D, _) Spec.Eval.Inputs.t -> Spec.Hyper.t array (*! EGPT[2] *)

        (** Optimization with the GNU Scientific library (GSL) *)
        module Gsl : sig
          (** [Optim_exception exn] is raised when an internal exception occurs,
              e.g. because GSL fails, or because a callback raised it. *)
          exception Optim_exception of exn

          (** [train ?step ?tol ?epsabs ?report_trained_model
              ?report_gradient_norm ?kernel ?sigma2 ?inducing ?n_rand_inducing
              ?learn_sigma2 ?hypers ~inputs ~targets ()] takes the optional
              initial optimizer step size [step], the optimizer line search
              tolerance [tol], the minimum gradient norm [epsabs] to achieve by
              the optimizer, callbacks for reporting intermediate results
              [report_trained_model] and [report_gradient_norm], an optional
              [kernel], noise level [sigma2], inducing inputs [inducing], number
              of randomly chosen inducing inputs [n_rand_inducing], a flag for
              whether the noise level should be learnt [learn_sigma2], an array
              of optional hyper parameters [hypers] which should be optimized,
              and the [inputs] and [targets].

              @return the trained model obtained by evidence maximization (=
              type II maximum likelihood).

              @param step default = [1e-1]
              @param tol default = [1e-1]
              @param epsabs default = [1e-1]
              @param report_trained_model default = ignored
              @param report_gradient_norm default = ignored
              @param kernel default = default kernel computed from specification
              @param sigma2 default = target variance
              @param inducing default = randomly selected subset of inputs
              @param n_rand_inducing default = 10% of inputs, at most 1000
              @param learn_sigma2 default = [true]
              @param hypers default = all hyper parameters
          *)
          val train :
            ?step : float ->
            ?tol : float ->
            ?epsabs : float ->
            ?report_trained_model : (iter : int -> ('D, 'd, 'm, 'n) Eval.Trained.t -> unit) -> (*! ITP *)
            ?report_gradient_norm : (iter : int -> float -> unit) ->
            kernel : ('D, 'd, 'm) Eval.Spec.Kernel.t -> (*! ITP *)
            ?sigma2 : float ->
            ?inducing : ('d, 'm) Eval.Spec.Inducing.t -> (*! ITP *)
            n_rand_inducing : 'm Slap.Size.t -> (*! ITP,O2L[1] *)
            ?learn_sigma2 : bool ->
            n_hypers : _ Slap.Size.t -> (*! EGPT[2] *)
            hypers : Spec.Hyper.t array -> (*! EGPT[2] *)
            inputs : ('D, 'n) Eval.Spec.Inputs.t -> (*! ITP *)
            targets : ('n, Slap.cnt) vec -> (*! ITP *)
            unit ->
            ('D, 'd, 'm, 'n) Eval.Trained.t (*! ITP *)
        end

        type ('D, 'd, 'm, 'n, 'ki, 'ko) learn_sigma2 (*! FT[1] *)
        val learn_sigma2 : ('D, 'd, 'm, 'n, 'k, 'k Slap.Size.s) learn_sigma2 (*! FT[1] *)
        val not_learn_sigma2 : ('D, 'd, 'm, 'n, 'k, 'k) learn_sigma2 (*! FT[1] *)

        module SGD : sig
          type ('D, 'd, 'm, 'n, 'ki, 'ko) t (*! ITP *)

          val create :
            ?tau : float ->
            ?eta0 : float ->
            ?step : int ->
            kernel : ('D, 'd, 'm) Eval.Spec.Kernel.t -> (*! ITP,O2L[2] *)
            ?sigma2 : float ->
            ?inducing : ('d, 'm) Eval.Spec.Inducing.t -> (*! ITP *)
            n_rand_inducing : 'm Slap.Size.t -> (*! ITP,O2L[1] *)
            learn_sigma2 : ('D, 'd, 'm, 'n, 'ki, 'ko) learn_sigma2 -> (*! FT[1] *)
            n_hypers : 'ki Slap.Size.t -> (*! EGPT[2] *)
            hypers : Spec.Hyper.t array -> (*! EGPT[2] *)
            inputs : ('D, 'n) Eval.Spec.Inputs.t -> (*! ITP *)
            targets : ('n, Slap.cnt) vec -> (*! ITP *)
            unit ->
            ('D, 'd, 'm, 'n, 'ki, 'ko) t (*! ITP *)

          val step : ('D, 'd, 'm, 'n, 'ki, 'ko) t -> ('D, 'd, 'm, 'n, 'ki, 'ko) t (*! ITP *)
          val gradient_norm : ('D, 'd, 'm, 'n, 'ki, 'ko) t -> float (*! ITP *)
          val get_trained : ('D, 'd, 'm, 'n, 'ki, 'ko) t -> ('D, 'd, 'm, 'n) Eval.Trained.t (*! ITP *)

          val get_eta : ('D, 'd, 'm, 'n, 'ki, 'ko) t -> float (*! ITP *)
          val get_step : ('D, 'd, 'm, 'n, 'ki, 'ko) t -> int (*! ITP *)

          val test :
            ?epsabs : float ->
            ?max_iter : int ->
            ?report :  (('D, 'd, 'm, 'n, 'ki, 'ko) t -> unit) -> (*! ITP *)
            ('D, 'd, 'm, 'n, 'ki, 'ko) t -> (*! ITP *)
            ('D, 'd, 'm, 'n, 'ki, 'ko) t (*! ITP *)
        end

        module SMD : sig
          type ('D, 'd, 'm, 'n, 'ki, 'ko) t (*! ITP *)

          val create :
            ?eps : float ->
            ?lambda : float ->
            ?mu : float ->
            ?eta0 : ('ko, Slap.cnt) vec -> (*! ITP *)
            ?nu0 : ('ko, Slap.cnt) vec -> (*! ITP *)
            kernel : ('D, 'd, 'm) Eval.Spec.Kernel.t -> (*! ITP,O2L[2] *)
            ?sigma2 : float ->
            ?inducing : ('d, 'm) Eval.Spec.Inducing.t -> (*! ITP *)
            n_rand_inducing : 'm Slap.Size.t -> (*! ITP,O2L[1] *)
            learn_sigma2 : ('D, 'd, 'm, 'n, 'ki, 'ko) learn_sigma2 -> (*! FT[1] *)
            n_hypers : 'ki Slap.Size.t -> (*! EGPT[2] *)
            hypers : Spec.Hyper.t array -> (*! EGPT[2] *)
            inputs : ('D, 'n) Eval.Spec.Inputs.t -> (*! ITP *)
            targets : ('n, Slap.cnt) vec -> (*! ITP *)
            unit ->
            ('D, 'd, 'm, 'n, 'ki, 'ko) t (*! ITP *)

          val step : ('D, 'd, 'm, 'n, 'ki, 'ko) t -> ('D, 'd, 'm, 'n, 'ki, 'ko) t (*! ITP *)
          val gradient_norm : ('D, 'd, 'm, 'n, 'ki, 'ko) t -> float (*! ITP *)
          val get_trained : ('D, 'd, 'm, 'n, 'ki, 'ko) t -> ('D, 'd, 'm, 'n) Eval.Trained.t (*! ITP *)

          val get_eta : ('D, 'd, 'm, 'n, 'ki, 'ko) t -> ('ko, Slap.cnt) vec (*! ITP *)
          val get_nu : ('D, 'd, 'm, 'n, 'ki, 'ko) t -> ('ko, Slap.cnt) vec (*! ITP *)

          val test :
            ?epsabs : float ->
            ?max_iter : int ->
            ?report : (('D, 'd, 'm, 'n, 'ki, 'ko) t -> unit) -> (*! ITP *)
            ('D, 'd, 'm, 'n, 'ki, 'ko) t -> (*! ITP *)
            ('D, 'd, 'm, 'n, 'ki, 'ko) t (*! ITP *)
        end
      end

(*
      (** Online learning *)
      module Online : sig
        type t

        val sgd :
          ?capacity : int ->
          ?eta0 : float -> ?tau : float -> Spec.Eval.Kernel.t -> t

        val smd :
          ?capacity : int ->
          ?eta0 : vec -> ?mu : float -> ?lam : float -> Spec.Eval.Kernel.t -> t

        val train : t -> Spec.Eval.Input.t -> target : float -> t

        val calc_mean_predictor : t -> Eval.Mean_predictor.t
        val calc_co_variance_predictor : t -> Eval.Co_variance_predictor.t
      end
*)

    end

  end

  (** Modules for global optimization *)
  module type Optimizer = sig

    (** Sub-modules for learning without derivatives. *)
    module Eval : Eval

    (** Sub-modules for global optimization. *)
    module Optimizer : sig
      module Spec : Specs.Optimizer with module Eval = Eval.Spec

      type t

      val create : ?max_memory : int -> (_, _,_) Spec.Eval.Kernel.t -> t (*! ITP *)

      val learn : t -> (_ Spec.Eval.Input.t * float) array -> t (*! ITP *)

      val calc_mpi_criterion : t -> _ Spec.Eval.Input.t -> float (*! ITP *)

      val calc_mpi_deriv : t -> _ Spec.Eval.Input.t (*! ITP *)
    end
  end
end
