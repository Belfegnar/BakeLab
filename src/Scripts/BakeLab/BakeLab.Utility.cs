#pragma warning disable
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Math = System.Math;

namespace Belfegnar.BakeLab {

	public struct Vector2i : System.IEquatable<Vector2i> {
		public int x;
		public int y;

		public Vector2i (int x, int y) {
			this.x = x;
			this.y = y;
		}

		/// <summary>
		/// Set x and y components of an existing Vector.
		/// </summary>
		public void Set (int x, int y) {
			this.x = x;
			this.y = y;
		}

		/// <summary>
		/// Access the /x/ or /y/ component using [0] or [1] respectively.
		/// </summary>
		public int this [int index] {
			get {
				switch (index) {
					case 0:
						return x;
					case 1:
						return y;
					default:
						throw new System.IndexOutOfRangeException (string.Format ("Invalid Vector2i index addressed: {0}!", index));
				}
			}

			set {
				switch (index) {
					case 0:
						x = value;
						break;
					case 1:
						y = value;
						break;
					default:
						throw new System.IndexOutOfRangeException (string.Format ("Invalid Vector2i index addressed: {0}!", index));
				}
			}
		}

		/// <summary>
		/// Returns the length of this vector (RO).
		/// </summary>
		public float Magnitude { get { return Mathf.Sqrt ((float) (x * x + y * y)); } }

		/// <summary>
		/// Returns the squared length of this vector (RO).
		/// </summary>
		public int SqrMagnitude { get { return x * x + y * y; } }

		/// <summary>
		/// Returns the distance between /a/ and /b/.
		/// </summary>
		public static float Distance (Vector2i a, Vector2i b) { return (a - b).Magnitude; }

		/// <summary>
		/// Returns a vector that is made from the smallest components of two vectors.
		/// </summary>
		public static Vector2i Min (Vector2i lhs, Vector2i rhs) { return new Vector2i (Mathf.Min (lhs.x, rhs.x), Mathf.Min (lhs.y, rhs.y)); }

		/// <summary>
		/// Returns a vector that is made from the largest components of two vectors.
		/// </summary>
		public static Vector2i Max (Vector2i lhs, Vector2i rhs) { return new Vector2i (Mathf.Max (lhs.x, rhs.x), Mathf.Max (lhs.y, rhs.y)); }

		/// <summary>
		/// Multiplies two vectors component-wise.
		/// </summary>
		public static Vector2i Scale (Vector2i a, Vector2i b) { return new Vector2i (a.x * b.x, a.y * b.y); }

		/// <summary>
		/// Multiplies every component of this vector by the same component of /scale/.
		/// </summary>
		public void Scale (Vector2i scale) { x *= scale.x; y *= scale.y; }

		public void Clamp (Vector2i min, Vector2i max) {
			x = min.x > x ? min.x : x;
			x = max.x < x ? max.x : x;
			y = min.y > y ? min.y : y;
			y = max.y < y ? max.y : y;
		}

		/// <summary>
		/// Converts a Vector2i to a [[Vector2]].
		/// </summary>
		public static implicit operator Vector2 (Vector2i v) {
			return new Vector2 (v.x, v.y);
		}

		/// <summary>
		/// Converts a Vector2i to a [[Vector3i]].
		/// </summary>
		public static explicit operator Vector3i (Vector2i v) {
			return new Vector3i (v.x, v.y, 0);
		}

		public static Vector2i FloorToInt (Vector2 v) {
			return new Vector2i (
				Mathf.FloorToInt (v.x),
				Mathf.FloorToInt (v.y)
			);
		}

		public static Vector2i CeilToInt (Vector2 v) {
			return new Vector2i (
				Mathf.CeilToInt (v.x),
				Mathf.CeilToInt (v.y)
			);
		}

		public static Vector2i RoundToInt (Vector2 v) {
			return new Vector2i (
				Mathf.RoundToInt (v.x),
				Mathf.RoundToInt (v.y)
			);
		}

		public static Vector2i operator + (Vector2i a, Vector2i b) {
			return new Vector2i (a.x + b.x, a.y + b.y);
		}

		public static Vector2i operator - (Vector2i a, Vector2i b) {
			return new Vector2i (a.x - b.x, a.y - b.y);
		}

		public static Vector2i operator * (Vector2i a, Vector2i b) {
			return new Vector2i (a.x * b.x, a.y * b.y);
		}

		public static Vector2i operator * (Vector2i a, int b) {
			return new Vector2i (a.x * b, a.y * b);
		}

		public static bool operator == (Vector2i lhs, Vector2i rhs) {
			return lhs.x == rhs.x && lhs.y == rhs.y;
		}

		public static bool operator != (Vector2i lhs, Vector2i rhs) {
			return !(lhs == rhs);
		}

		public override bool Equals (object other) {
			if (!(other is Vector2i)) return false;

			return Equals ((Vector2i) other);
		}

		public bool Equals (Vector2i other) {
			return x.Equals (other.x) && y.Equals (other.y);
		}

		public override int GetHashCode () {
			return x.GetHashCode () ^ (y.GetHashCode () << 2);
		}

		public override string ToString () {
			return string.Format ("({0}, {1})", x, y);
		}

		private static readonly Vector2i zero = new Vector2i (0, 0);
		private static readonly Vector2i one = new Vector2i (1, 1);
		private static readonly Vector2i up = new Vector2i (0, 1);
		private static readonly Vector2i down = new Vector2i (0, -1);
		private static readonly Vector2i left = new Vector2i (-1, 0);
		private static readonly Vector2i right = new Vector2i (1, 0);
	}

	public struct Vector3i : System.IEquatable<Vector3i> {
		public int x;
		public int y;
		public int z;

		public Vector3i (int x, int y, int z) {
			this.x = x;
			this.y = y;
			this.z = z;
		}

		/// <summary>
		/// Set x, y and z components of an existing Vector.
		/// </summary>
		public void Set (int x, int y, int z) {
			this.x = x;
			this.y = y;
			this.z = z;
		}

		/// <summary>
		/// Access the /x/, /y/ or /z/ component using [0], [1] or [2] respectively.
		/// </summary>
		public int this [int index] {
			get {
				switch (index) {
					case 0:
						return x;
					case 1:
						return y;
					case 2:
						return z;
					default:
						throw new System.IndexOutOfRangeException (string.Format ("Invalid Vector3i index addressed: {0}!", index));
				}
			}

			set {
				switch (index) {
					case 0:
						x = value;
						break;
					case 1:
						y = value;
						break;
					case 2:
						z = value;
						break;
					default:
						throw new System.IndexOutOfRangeException (string.Format ("Invalid Vector3i index addressed: {0}!", index));
				}
			}
		}

		/// <summary>
		/// Returns the length of this vector (RO).
		/// </summary>
		public float Magnitude { get { return Mathf.Sqrt ((float) (x * x + y * y + z * z)); } }

		/// <summary>
		/// Returns the squared length of this vector (RO).
		/// </summary>
		public int SqrMagnitude { get { return x * x + y * y + z * z; } }

		/// <summary>
		/// Returns the distance between /a/ and /b/.
		/// </summary>
		public static float Distance (Vector3i a, Vector3i b) { return (a - b).Magnitude; }

		/// <summary>
		/// Returns a vector that is made from the smallest components of two vectors.
		/// </summary>
		public static Vector3i Min (Vector3i lhs, Vector3i rhs) { return new Vector3i (Mathf.Min (lhs.x, rhs.x), Mathf.Min (lhs.y, rhs.y), Mathf.Min (lhs.z, rhs.z)); }

		/// <summary>
		/// Returns a vector that is made from the largest components of two vectors.
		/// </summary>
		public static Vector3i Max (Vector3i lhs, Vector3i rhs) { return new Vector3i (Mathf.Max (lhs.x, rhs.x), Mathf.Max (lhs.y, rhs.y), Mathf.Max (lhs.z, rhs.z)); }

		/// <summary>
		/// Multiplies two vectors component-wise.
		/// </summary>
		public static Vector3i Scale (Vector3i a, Vector3i b) { return new Vector3i (a.x * b.x, a.y * b.y, a.z * b.z); }

		/// <summary>
		/// Multiplies every component of this vector by the same component of /scale/.
		/// </summary>
		public void Scale (Vector3i scale) { x *= scale.x; y *= scale.y; z *= scale.z; }

		public void Clamp (Vector3i min, Vector3i max) {
			x = Mathf.Max (min.x, x);
			x = Mathf.Min (max.x, x);
			y = Mathf.Max (min.y, y);
			y = Mathf.Min (max.y, y);
			z = Mathf.Max (min.z, z);
			z = Mathf.Min (max.z, z);
		}

		/// <summary>
		/// Converts a Vector3i to a [[Vector3]].
		/// </summary>
		public static implicit operator Vector3 (Vector3i v) {
			return new Vector3 (v.x, v.y, v.z);
		}

		/// <summary>
		/// Converts a Vector3i to a [[Vector2i]].
		/// </summary>
		public static explicit operator Vector2i (Vector3i v) {
			return new Vector2i (v.x, v.y);
		}

		public static Vector3i FloorToInt (Vector3 v) {
			return new Vector3i (
				Mathf.FloorToInt (v.x),
				Mathf.FloorToInt (v.y),
				Mathf.FloorToInt (v.z)
			);
		}

		public static Vector3i CeilToInt (Vector3 v) {
			return new Vector3i (
				Mathf.CeilToInt (v.x),
				Mathf.CeilToInt (v.y),
				Mathf.CeilToInt (v.z)
			);
		}

		public static Vector3i RoundToInt (Vector3 v) {
			return new Vector3i (
				Mathf.RoundToInt (v.x),
				Mathf.RoundToInt (v.y),
				Mathf.RoundToInt (v.z)
			);
		}

		public static Vector3i operator + (Vector3i a, Vector3i b) {
			return new Vector3i (a.x + b.x, a.y + b.y, a.z + b.z);
		}

		public static Vector3i operator - (Vector3i a, Vector3i b) {
			return new Vector3i (a.x - b.x, a.y - b.y, a.z - b.z);
		}

		public static Vector3i operator * (Vector3i a, Vector3i b) {
			return new Vector3i (a.x * b.x, a.y * b.y, a.z * b.z);
		}

		public static Vector3i operator * (Vector3i a, int b) {
			return new Vector3i (a.x * b, a.y * b, a.z * b);
		}

		public static bool operator == (Vector3i lhs, Vector3i rhs) {
			return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
		}

		public static bool operator != (Vector3i lhs, Vector3i rhs) {
			return !(lhs == rhs);
		}

		public override bool Equals (object other) {
			if (!(other is Vector3i)) return false;

			return Equals ((Vector3i) other);
		}

		public bool Equals (Vector3i other) {
			return this == other;
		}

		public override int GetHashCode () {
			var yHash = y.GetHashCode ();
			var zHash = z.GetHashCode ();
			return x.GetHashCode () ^ (yHash << 4) ^ (yHash >> 28) ^ (zHash >> 4) ^ (zHash << 28);
		}

		public override string ToString () {
			return string.Format ("({0}, {1}, {2})", x, y, z);
		}

		public string ToString (string format) {
			return string.Format ("({0}, {1}, {2})", x.ToString (format), y.ToString (format), z.ToString (format));
		}

		public static Vector3i zero { get { return s_Zero; } }
		public static Vector3i one { get { return s_One; } }
		public static Vector3i up { get { return s_Up; } }
		public static Vector3i down { get { return s_Down; } }
		public static Vector3i left { get { return s_Left; } }
		public static Vector3i right { get { return s_Right; } }

		private static readonly Vector3i s_Zero = new Vector3i (0, 0, 0);
		private static readonly Vector3i s_One = new Vector3i (1, 1, 1);
		private static readonly Vector3i s_Up = new Vector3i (0, 1, 0);
		private static readonly Vector3i s_Down = new Vector3i (0, -1, 0);
		private static readonly Vector3i s_Left = new Vector3i (-1, 0, 0);
		private static readonly Vector3i s_Right = new Vector3i (1, 0, 0);
	}

	public struct Vector2d {

		/// <summary>
		/// x coordinate.
		/// </summary>
		public double x;

		/// <summary>
		/// x coordinate.
		/// </summary>
		public double y;

		public Vector2d (double xd, double yd) {
			x = xd;
			y = yd;
		}

		public static bool operator == (Vector2d lhs, Vector2d v) {
			return (lhs.x == v.x && lhs.y == v.y);
		}

		public static bool operator != (Vector2d lhs, Vector2d v) {
			return (lhs.x != v.x || lhs.y != v.y);
		}

		/// <summary>
		/// Is instance equals with specified one.
		/// </summary>
		/// <param name="rhs">Specified instance for comparation.</param>
		public override bool Equals (object rhs) {
			if (!(rhs is Vector2d)) {
				return false;
			}
			return this == (Vector2d) rhs;
		}

		/// <summary>
		/// Return formatted X/Y/Z values.
		/// </summary>
		public override string ToString () {
			return string.Format ("({0}, {1})", x, y);
		}

		/// <summary>
		/// Vectors hash code. 
		/// </summary>
		public override int GetHashCode () {
			double hashcode = 23;
			hashcode = (hashcode * 37) + x;
			hashcode = (hashcode * 37) + y;

			return (int) hashcode;
		}

		public static Vector2d operator + (Vector2d lhs, Vector2d v) {
			return new Vector2d (lhs.x + v.x, lhs.y + v.y);
		}

		public static Vector2d operator - (Vector2d lhs, Vector2d v) {
			return new Vector2d (lhs.x - v.x, lhs.y - v.y);
		}

		public static Vector2d operator * (Vector2d lhs, Vector2d v) {
			return new Vector2d (lhs.x * v.x, lhs.y * v.y);
		}

		public static Vector2d operator * (Vector2d lhs, double scalar) {
			return new Vector2d (lhs.x * scalar, lhs.y * scalar);
		}

		public static Vector2d operator / (Vector2d lhs, Vector2d v) {
			return new Vector2d (lhs.x / v.x, lhs.y / v.y);
		}

		public static Vector2d operator / (Vector2d lhs, double scalar) {
#if DEBUG
			Debug.Assert (scalar != 0.0);
#endif
			double inv = 1.0 / scalar;
			return new Vector2d (lhs.x * inv, lhs.y * inv);
		}

		public double Length () {
			return System.Math.Sqrt (x * x + y * y);
		}

		public Vector2d Normalize () {
			double length = System.Math.Sqrt (x * x + y * y);
			double invLength = 1.0 / length;
			return new Vector2d (x * invLength, y * invLength);
		}

		public Vector2 CastVector2 () {
			return new Vector2 ((float) x, (float) y);
		}

		public double Magnitude {
			get {
				return System.Math.Sqrt (x * x + y * y);
			}
		}

		public double SqrMagnitude {
			get {
				return x * x + y * y;
			}
		}

		public static readonly Vector2d ZERO = new Vector2d (0, 0);

		public static readonly Vector2d UNIT_X = new Vector2d (1, 0);

		public static readonly Vector2d UNIT_Y = new Vector2d (0, 1);
	}

	public struct Vector3d {
		/// <summary>
		/// x coordinate.
		/// </summary>
		public double x;

		/// <summary>
		/// x coordinate.
		/// </summary>
		public double y;

		/// <summary>
		/// z coordinate.
		/// </summary>
		public double z;

		public Vector3d (double xd, double yd, double zd) {
			x = xd;
			y = yd;
			z = zd;
		}

		public static implicit operator Vector3d (double[] op) {
			Vector3d result;
			result.x = op[0];
			result.y = op[1];
			result.z = op[2];
			return result;
		}

		public double this [int index] {
			get {
				if (index == 0) {
					return x;
				} else if (index == 1) {
					return y;
				} else if (index == 2) {
					return z;
				}
				throw new System.IndexOutOfRangeException ();
				//			return 0;
			}
			set {
				if (index == 0) {
					x = value;
					return;
				} else if (index == 1) {
					y = value;
					return;
				} else if (index == 2) {
					z = value;
					return;
				}
				throw new System.IndexOutOfRangeException ();
			}
		}

		public static bool operator == (Vector3d lhs, Vector3d v) {
			return (lhs.x == v.x && lhs.y == v.y && lhs.z == v.z);
		}

		public static bool operator != (Vector3d lhs, Vector3d v) {
			return (lhs.x != v.x || lhs.y != v.y || lhs.z != v.z);
		}

		/// <summary>
		/// Is instance equals with specified one.
		/// </summary>
		/// <param name="rhs">Specified instance for comparation.</param>
		public override bool Equals (object rhs) {
			if (!(rhs is Vector3d)) {
				return false;
			}
			return this == (Vector3d) rhs;
		}

		/// <summary>
		/// Return formatted X/Y/Z values.
		/// </summary>
		public override string ToString () {
			return string.Format ("({0}, {1}, {2})", x, y, z);
		}

		/// <summary>
		/// Vectors hash code. 
		/// </summary>
		public override int GetHashCode () {
			double hashcode = 23;
			hashcode = (hashcode * 37) + x;
			hashcode = (hashcode * 37) + y;
			hashcode = (hashcode * 37) + z;

			return (int) hashcode;
		}

		public static Vector3d operator + (Vector3d v1, Vector3d v2) {
			Vector3d result;
			result.x = v1.x + v2.x;
			result.y = v1.y + v2.y;
			result.z = v1.z + v2.z;
			return result;
		}

		public static Vector3d operator - (Vector3d v1, Vector3d v2) {
			Vector3d result;
			result.x = v1.x - v2.x;
			result.y = v1.y - v2.y;
			result.z = v1.z - v2.z;
			return result;
		}

		public static Vector3d operator * (Vector3d lhs, Vector3d v) {
			return new Vector3d (lhs.x * v.x, lhs.y * v.y, lhs.z * v.z);
		}

		public static Vector3d operator * (Vector3d v1, double scalar) {
			Vector3d result;
			result.x = v1.x * scalar;
			result.y = v1.y * scalar;
			result.z = v1.z * scalar;
			return result;
		}

		public static Vector3d operator * (double scalar, Vector3d v1) {
			Vector3d result;
			result.x = v1.x * scalar;
			result.y = v1.y * scalar;
			result.z = v1.z * scalar;
			return result;
		}

		public static Vector3d operator / (Vector3d lhs, Vector3d v) {
			return new Vector3d (lhs.x / v.x, lhs.y / v.y, lhs.z / v.z);
		}

		public static Vector3d operator / (Vector3d v1, double scalar) {
#if DEBUG
			Debug.Assert (scalar != 0.0);
#endif
			double inv = 1.0 / scalar;
			Vector3d result;
			result.x = v1.x * inv;
			result.y = v1.y * inv;
			result.z = v1.z * inv;
			return result;
		}

		public double Length () {
			return System.Math.Sqrt (x * x + y * y + z * z);
		}

		public double SquaredLength () {
			return x * x + y * y + z * z;
		}

		public double DistanceTo (Vector3d v) {
			double xx = x - v.x;
			double yy = y - v.y;
			double zz = z - v.z;
			return System.Math.Sqrt (xx * xx + yy * yy + zz * zz);
		}

		public double SquaredDistanceTo (Vector3d v) {
			double xx = x - v.x;
			double yy = y - v.y;
			double zz = z - v.z;
			return xx * xx + yy * yy + zz * zz;
		}

		public double Dotproduct (Vector3d v) {
			return (x * v.x + y * v.y + z * v.z);
		}

		public Vector3d Normalize () {
			double length = System.Math.Sqrt (x * x + y * y + z * z);
			double inv = 1.0 / length;
			Vector3d result;
			result.x = x * inv;
			result.y = y * inv;
			result.z = z * inv;
			return result;
		}

		public Vector3d Normalize (double l) {
			double length = Math.Sqrt (x * x + y * y + z * z);
			double inv = l / length;
			Vector3d result;
			result.x = x * inv;
			result.y = y * inv;
			result.z = z * inv;
			return result;
		}

		public Vector3d Normalize (out double previousLength) {
			previousLength = Math.Sqrt (x * x + y * y + z * z);
			double inv = 1.0 / previousLength;
			Vector3d result;
			result.x = x * inv;
			result.y = y * inv;
			result.z = z * inv;
			return result;
		}

		public Vector3d Normalize (double l, out double previousLength) {
			previousLength = Math.Sqrt (x * x + y * y + z * z);
			double inv = l / previousLength;
			Vector3d result;
			result.x = x * inv;
			result.y = y * inv;
			result.z = z * inv;
			return result;
		}

		public Vector3d CrossProduct (Vector3d v) {
			Vector3d result;
			result.x = y * v.z - z * v.y;
			result.y = z * v.x - x * v.z;
			result.z = x * v.y - y * v.x;
			return result;
		}

		public static Vector3d operator % (Vector3d u, Vector3d v) {
			Vector3d result;
			result.x = u.y * v.z - u.z * v.y;
			result.y = u.z * v.x - u.x * v.z;
			result.z = u.x * v.y - u.y * v.x;
			return result;
		}

		public Vector2d xy () {
			Vector2d v;
			v.x = x;
			v.y = y;
			return v;
		}

		public Vector3 ToVector3 () {
			Vector3 result;
			result.x = (float) x;
			result.y = (float) y;
			result.z = (float) z;
			return result;
		}

		public static readonly Vector3d ZERO = new Vector3d (0, 0, 0);

		public static readonly Vector3d UNIT_X = new Vector3d (1, 0, 0);

		public static readonly Vector3d UNIT_Y = new Vector3d (0, 1, 0);

		public static readonly Vector3d UNIT_Z = new Vector3d (0, 0, 1);
	}
}

namespace MathNet.Numerics.LinearAlgebra.Double.Factorization {
	using System;
	using CSparse.Double;
	using CSparse;
	using MathNet.Numerics.LinearAlgebra.Storage;
	using MathNet.Numerics.Properties;

	// Create an alias for CSparse's SparseLU class.
	using CSparseLU = CSparse.Double.Factorization.SparseLU;
	using CSparseCholesky = CSparse.Double.Factorization.SparseCholesky;
	using CSparseLDL = CSparse.Double.Factorization.SparseLDL;

	public class SparseLU {
		int n;
		CSparseLU lu;

		private SparseLU (CSparseLU lu, int n) {
			this.n = n;
			this.lu = lu;
		}

		/// <summary>
		/// Compute the sparse LU factorization for given matrix.
		/// </summary>
		/// <param name="matrix">The matrix to factorize.</param>
		/// <param name="ordering">The column ordering method to use.</param>
		/// <param name="tol">Partial pivoting tolerance (form 0.0 to 1.0).</param>
		/// <returns>Sparse LU factorization.</returns>
		public static SparseLU Create (MathNet.Numerics.LinearAlgebra.Double.SparseMatrix matrix, ColumnOrdering ordering,
			double tol = 1.0) {
			int n = matrix.RowCount;

			// Check for proper dimensions.
			if (n != matrix.ColumnCount) {
				throw new ArgumentException (Resources.MatrixMustBeSparse);
			}

			// Get CSR storage.
			var storage = (SparseCompressedRowMatrixStorage<double>) matrix.Storage;

			// Create CSparse matrix.
			var A = new CompressedColumnStorage (n, n);

			// Assign storage arrays.
			A.ColumnPointers = storage.RowPointers;
			A.RowIndices = storage.ColumnIndices;
			A.Values = storage.Values;

			return new SparseLU (CSparseLU.Create (A, ordering, tol), n);
		}

		/// <summary>
		/// Solves a system of linear equations, <c>Ax = b</c>, with A LU factorized.
		/// </summary>
		/// <param name="input">The right hand side vector, <c>b</c>.</param>
		/// <param name="result">The left hand side vector, <c>x</c>.</param>
		public void Solve (Vector<double> input, Vector<double> result) {
			// Check for proper arguments.
			if (input == null) {
				throw new ArgumentNullException ("input");
			}

			if (result == null) {
				throw new ArgumentNullException ("result");
			}

			// Check for proper dimensions.
			if (input.Count != result.Count) {
				throw new ArgumentException (Resources.ArgumentVectorsSameLength);
			}

			if (input.Count != n) {
				throw new ArgumentException ("Dimensions don't match", "input");
			}

			var b = input.Storage as DenseVectorStorage<double>;
			var x = result.Storage as DenseVectorStorage<double>;

			if (b == null || x == null) {
				throw new NotSupportedException ("Expected dense vector storage.");
			}

			lu.SolveTranspose (b.Data, x.Data);
		}
	}

	public class SparseCholesky {
		int n;
		CSparseCholesky lu;

		private SparseCholesky (CSparseCholesky lu, int n) {
			this.n = n;
			this.lu = lu;
		}

		/// <summary>
		/// Compute the sparse LU factorization for given matrix.
		/// </summary>
		/// <param name="matrix">The matrix to factorize.</param>
		/// <param name="ordering">The column ordering method to use.</param>
		/// <param name="tol">Partial pivoting tolerance (form 0.0 to 1.0).</param>
		/// <returns>Sparse LU factorization.</returns>
		public static SparseCholesky Create (MathNet.Numerics.LinearAlgebra.Double.SparseMatrix matrix, ColumnOrdering ordering,
			double tol = 1.0) {
			int n = matrix.RowCount;

			// Check for proper dimensions.
			if (n != matrix.ColumnCount) {
				throw new ArgumentException (Resources.MatrixMustBeSparse);
			}

			// Get CSR storage.
			var storage = (SparseCompressedRowMatrixStorage<double>) matrix.Storage;

			// Create CSparse matrix.
			var A = new CompressedColumnStorage (n, n);

			// Assign storage arrays.
			A.ColumnPointers = storage.RowPointers;
			A.RowIndices = storage.ColumnIndices;
			A.Values = storage.Values;

			return new SparseCholesky (CSparseCholesky.Create (A, ordering /*, tol*/ ), n);
		}

		/// <summary>
		/// Solves a system of linear equations, <c>Ax = b</c>, with A LU factorized.
		/// </summary>
		/// <param name="input">The right hand side vector, <c>b</c>.</param>
		/// <param name="result">The left hand side vector, <c>x</c>.</param>
		public void Solve (Vector<double> input, Vector<double> result) {
			// Check for proper arguments.
			if (input == null) {
				throw new ArgumentNullException ("input");
			}

			if (result == null) {
				throw new ArgumentNullException ("result");
			}

			// Check for proper dimensions.
			if (input.Count != result.Count) {
				throw new ArgumentException (Resources.ArgumentVectorsSameLength);
			}

			if (input.Count != n) {
				throw new ArgumentException ("Dimensions don't match", "input");
			}

			var b = input.Storage as DenseVectorStorage<double>;
			var x = result.Storage as DenseVectorStorage<double>;

			if (b == null || x == null) {
				throw new NotSupportedException ("Expected dense vector storage.");
			}

			lu.Solve (b.Data, x.Data);
		}

		public void Solve (double[] input, double[] result) {
			lu.Solve (input, result);
		}
	}

	public class SparseLDL {
		int n;
		CSparseLDL lu;

		private SparseLDL (CSparseLDL lu, int n) {
			this.n = n;
			this.lu = lu;
		}

		/// <summary>
		/// Compute the sparse LU factorization for given matrix.
		/// </summary>
		/// <param name="matrix">The matrix to factorize.</param>
		/// <param name="ordering">The column ordering method to use.</param>
		/// <param name="tol">Partial pivoting tolerance (form 0.0 to 1.0).</param>
		/// <returns>Sparse LU factorization.</returns>
		public static SparseLDL Create (MathNet.Numerics.LinearAlgebra.Double.SparseMatrix matrix, ColumnOrdering ordering,
			double tol = 1.0) {
			int n = matrix.RowCount;

			// Check for proper dimensions.
			if (n != matrix.ColumnCount) {
				throw new ArgumentException (Resources.MatrixMustBeSparse);
			}

			// Get CSR storage.
			var storage = (SparseCompressedRowMatrixStorage<double>) matrix.Storage;

			// Create CSparse matrix.
			var A = new CompressedColumnStorage (n, n);

			// Assign storage arrays.
			A.ColumnPointers = storage.RowPointers;
			A.RowIndices = storage.ColumnIndices;
			A.Values = storage.Values;

			return new SparseLDL (new CSparseLDL (A, ordering /*, tol*/ ), n);
		}

		/// <summary>
		/// Solves a system of linear equations, <c>Ax = b</c>, with A LU factorized.
		/// </summary>
		/// <param name="input">The right hand side vector, <c>b</c>.</param>
		/// <param name="result">The left hand side vector, <c>x</c>.</param>
		public void Solve (Vector<double> input, Vector<double> result) {
			// Check for proper arguments.
			if (input == null) {
				throw new ArgumentNullException ("input");
			}

			if (result == null) {
				throw new ArgumentNullException ("result");
			}

			// Check for proper dimensions.
			if (input.Count != result.Count) {
				throw new ArgumentException (Resources.ArgumentVectorsSameLength);
			}

			if (input.Count != n) {
				throw new ArgumentException ("Dimensions don't match", "input");
			}

			var b = input.Storage as DenseVectorStorage<double>;
			var x = result.Storage as DenseVectorStorage<double>;

			if (b == null || x == null) {
				throw new NotSupportedException ("Expected dense vector storage.");
			}

			lu.Solve (b.Data, x.Data);
		}

		public void Solve (double[] input, double[] result) {
			lu.Solve (input, result);
		}
	}
}

namespace MathNet.Numerics.LinearAlgebra.Double {
	using CSparse;

	using LU = MathNet.Numerics.LinearAlgebra.Double.Factorization.SparseLU;
	using Cholesky = MathNet.Numerics.LinearAlgebra.Double.Factorization.SparseCholesky;
	using LDL = MathNet.Numerics.LinearAlgebra.Double.Factorization.SparseLDL;

	public static class SparseMatrixExtensions {
		/// <summary>
		/// Compute the sparse LU factorization for given matrix.
		/// </summary>
		/// <param name="matrix">The matrix to factorize.</param>
		/// <param name="tol">Partial pivoting tolerance (form 0.0 to 1.0).</param>
		/// <returns></returns>
		public static LU SparseLU (this SparseMatrix matrix, double tol = 1.0) {
			return LU.Create (matrix, ColumnOrdering.MinimumDegreeAtPlusA, tol);
		}

		public static Cholesky SparseCholesky (this SparseMatrix matrix, double tol = 1.0) {
			return Cholesky.Create (matrix, ColumnOrdering.MinimumDegreeAtPlusA, tol);
		}

		public static LDL SparseLDL (this SparseMatrix matrix, double tol = 1.0) {
			return LDL.Create (matrix, ColumnOrdering.MinimumDegreeAtPlusA, tol);
		}
	}
}