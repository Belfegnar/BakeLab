using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using EdgeMap = System.Collections.Generic.Dictionary<BelfegnarInc.BakeLab.Vector2i, BelfegnarInc.BakeLab.Butterfly>;
using TripletMap = System.Collections.Generic.Dictionary<BelfegnarInc.BakeLab.Vector2i, System.Double>;
using Timer = System.Diagnostics.Stopwatch;
using MathNet.Numerics.LinearAlgebra.Double;

namespace BelfegnarInc.BakeLab {
	public static class Filter {

		const double eps = 0.0;

		public static void MapSignalToVertices (
			Scene scene,
			int[] numSamplesPerInstance,
			Samples samples,
			float[] values,
			VertexFilterMode mode,
			float regularizationWeight,
			float[][] vertexSignals,
			bool parallel = false
		) {
			if (mode == VertexFilterMode.AREA_BASED) {
				DoFilter (scene, numSamplesPerInstance, samples, values, vertexSignals, parallel);
			} else if (mode == VertexFilterMode.LEAST_SQUARES) {
				DoFilterLeastSquares (scene, numSamplesPerInstance, samples, values, regularizationWeight, vertexSignals, parallel);
			}
#if DEBUG
			else {
				Debug.Assert (false, "Invalid vertex filter mode!");
			}
#endif
		}

		#region area_based

		static void DoFilter (
			Scene scene,
			int[] numSamplesPerInstance,
			Samples samples,
			float[] values,
			float[][] vertexSignals,
			bool parallel = false
		) {
			int[] sampleOffsetPerInstance = new int[scene.numInstances]; {
				int sample_offset = 0;
				for (int i = 0; i < scene.numInstances; ++i) {
					sampleOffsetPerInstance[i] = sample_offset;
					sample_offset += numSamplesPerInstance[i];
				}
			}

			if (parallel) {
				//#pragma omp parallel for
				System.Threading.Tasks.Parallel.For (0, scene.numInstances, (i) => {
					int sampleOffset = sampleOffsetPerInstance[i];

					// Point to samples for this instance
					Samples instanceSamples = new Samples ();
					instanceSamples.numSamples = numSamplesPerInstance[i];
					instanceSamples.samplePositions = samples.samplePositions;
					instanceSamples.sampleNormals = samples.sampleNormals;
					instanceSamples.sampleFaceNormals = samples.sampleFaceNormals;
					instanceSamples.sampleInfos = samples.sampleInfos;
					instanceSamples.sampleOffset = sampleOffset;

					FilterMeshAreaWeighted (scene.meshes[scene.instances[i].meshIndex], instanceSamples, values, vertexSignals[i]);
				});
			} else {
				for (int i = 0; i < scene.numInstances; ++i) {
					int sampleOffset = sampleOffsetPerInstance[i];

					// Point to samples for this instance
					Samples instanceSamples = new Samples ();
					instanceSamples.numSamples = numSamplesPerInstance[i];
					instanceSamples.samplePositions = samples.samplePositions;
					instanceSamples.sampleNormals = samples.sampleNormals;
					instanceSamples.sampleFaceNormals = samples.sampleFaceNormals;
					instanceSamples.sampleInfos = samples.sampleInfos;
					instanceSamples.sampleOffset = sampleOffset;

					FilterMeshAreaWeighted (scene.meshes[scene.instances[i].meshIndex], instanceSamples, values, vertexSignals[i]);
				}
			}
		}

		/// <summary>
		/// Splat area-weighted samples onto vertices
		/// </summary>
		static void FilterMeshAreaWeighted (
			Mesh mesh,
			Samples samples,
			float[] values,
			float[] vertexSignals
		) {
			int sampleOffset = samples.sampleOffset;
			double[] weights = new double[mesh.numVertices];
			for (int i = 0; i < mesh.numVertices; i++)
				vertexSignals[i] = 0.0f;

			int[] triangles = mesh.triangles;
			Vector3i tri = new Vector3i ();
			for (int i = 0; i < samples.numSamples; ++i) {
				SampleInfo info = samples.sampleInfos[i + sampleOffset];
				tri.x = triangles[info.triIdx * 3];
				tri.y = triangles[info.triIdx * 3 + 1];
				tri.z = triangles[info.triIdx * 3 + 2];

				float val = values[i + sampleOffset];
				Vector3 w = new Vector3 (info.dA * info.bary.x,
					info.dA * info.bary.y,
					info.dA * info.bary.z);

				vertexSignals[tri.x] += w.x * val;
				weights[tri.x] += w.x;
				vertexSignals[tri.y] += w.y * val;
				weights[tri.y] += w.y;
				vertexSignals[tri.z] += w.z * val;
				weights[tri.z] += w.z;
			}

			// Normalize
			for (int k = 0; k < mesh.numVertices; ++k) {
				if (weights[k] > 0.0) vertexSignals[k] /= (float) (weights[k]);
			}
		}

		#endregion

		#region least_squares

		static double TriangleArea (Vector3d a, Vector3d b, Vector3d c) {
			Vector3d ba = b - a, ca = c - a;
			Vector3d crop = ba.crossProduct (ca);
			return crop.length () * 0.5;
		}

		/// <summary>
		/// Embeds 3D triangle v[0], v[1], v[2] into a plane, such that:
		///  p[0] = (0, 0), p[1] = (0, positive number), p[2] = (positive number, any number)
		/// If triangle is close to degenerate returns false and p is undefined.
		/// </summary>
		static bool PlanarizeTriangle (Vector3d[] v, Vector2d[] p) {
			double l01 = (v[0] - v[1]).length ();
			double l02 = (v[0] - v[2]).length ();
			double l12 = (v[1] - v[2]).length ();

			if (l01 <= eps || l02 <= eps || l12 <= eps) return false;

			double p2y = (l02 * l02 + l01 * l01 - l12 * l12) / (2.0 * l01);
			double tmp1 = l02 * l02 - p2y * p2y;
			if (tmp1 <= eps) return false;

			p[0] = new Vector2d (0.0, 0.0);
			p[1] = new Vector2d (0.0, l01);
			p[2] = new Vector2d (System.Math.Sqrt (tmp1), p2y);
			return true;
		}

		/// <summary>
		/// Computes gradient operator (2 x 3 matrix 'grad') for a planar triangle.  If
		/// 'normalize' is false then division by determinant is off (and thus the
		/// routine cannot fail even for degenerate triangles).
		/// </summary>
		static bool TriGrad2D (Vector2d[] p, bool normalize, ref /*Matrix23*/ MathNet.Numerics.LinearAlgebra.Matrix<double> grad) {
			double det = 1.0;
			if (normalize) {
				det = -p[0].y * p[1].x + p[0].x * p[1].y + p[0].y * p[2].x -
					p[1].y * p[2].x - p[0].x * p[2].y + p[1].x * p[2].y;

				if (System.Math.Abs (det) <= eps) {
					return false;
				}
			}

			grad[0, 0] = p[1].y - p[2].y;
			grad[0, 1] = p[2].y - p[0].y;
			grad[0, 2] = p[0].y - p[1].y;

			grad[1, 0] = p[2].x - p[1].x;
			grad[1, 1] = p[0].x - p[2].x;
			grad[1, 2] = p[1].x - p[0].x;

			grad.Multiply (1f / det, grad);
			//grad /= det;
			return true;
		}

		/// <summary>
		/// Computes difference of gradients operator (2 x 4 matrix 'GD') for a butterfly, i.e., 
		/// two edge-adjacent triangles.
		/// Performs normalization so that units are [m], so GD^T * GD will have units of area [m^2]:
		/// </summary>
		static bool ButterflyGradDiff (Vector3d[] v, ref /*Matrix24*/ MathNet.Numerics.LinearAlgebra.Matrix<double> GD) {
			Vector3d[] v1 = new Vector3d[] { v[0], v[1], v[2] };
			Vector3d[] v2 = new Vector3d[] { v[0], v[1], v[3] };
			Vector2d[] p1 = new Vector2d[3], p2 = new Vector2d[3];
			bool success = PlanarizeTriangle (v1, p1);
			if (!success) return false;
			success = PlanarizeTriangle (v2, p2);
			if (!success) return false;
			p2[2].x *= -1.0; // flip the x coordinate of the last vertex of the second triangle so we get a butterfly

			/*Matrix23*/
			MathNet.Numerics.LinearAlgebra.Matrix<double> grad1 = MathNet.Numerics.LinearAlgebra.Double.Matrix.Build.Dense (2, 3), grad2 = MathNet.Numerics.LinearAlgebra.Double.Matrix.Build.Dense (2, 3);
			success = TriGrad2D (p1, /*normalize*/ true, ref grad1);
			if (!success) return false;
			success = TriGrad2D (p2, true, ref grad2);
			if (!success) return false;

			/*Matrix24*/
			MathNet.Numerics.LinearAlgebra.Matrix<double> gradExt1 = MathNet.Numerics.LinearAlgebra.Double.Matrix.Build.Dense (2, 4), gradExt2 = MathNet.Numerics.LinearAlgebra.Double.Matrix.Build.Dense (2, 4);
			for (int i = 0; i < 2; i++) {
				gradExt1[i, 0] = grad1[i, 0];
				gradExt1[i, 1] = grad1[i, 1];
				gradExt1[i, 2] = grad1[i, 2];
				gradExt1[i, 3] = 0.0;
				gradExt2[i, 0] = grad2[i, 0];
				gradExt2[i, 1] = grad2[i, 1];
				gradExt2[i, 2] = 0.0;
				gradExt2[i, 3] = grad2[i, 2];
			}
			GD = gradExt1 - gradExt2;

			double area1 = TriangleArea (v1[0], v1[1], v1[2]);
			double area2 = TriangleArea (v2[0], v2[1], v2[2]);
			GD *= (area1 + area2);

			return true;
		}

		static void BuildRegularizationMatrix (
			Mesh mesh,
			ref MathNet.Numerics.LinearAlgebra.Double.SparseMatrix regularizationMatrix,
			Timer timer
		) {
			timer.Start ();
			if (mesh.regularizationMatrix != null) {
				regularizationMatrix = mesh.regularizationMatrix;
			} else {
				EdgeBasedRegularizer (mesh.vertices, mesh.numVertices, mesh.triangles, mesh.numTriangles, ref regularizationMatrix);
				mesh.regularizationMatrix = regularizationMatrix;
			}
			timer.Stop ();
		}

		static void EdgeBasedRegularizer (
			Vector3[] verts,
			int numVerts,
			int[] faces,
			int numFaces,
			ref MathNet.Numerics.LinearAlgebra.Double.SparseMatrix regularizationMatrix) {

			// Build edge map.  Each non-boundary edge stores the two opposite "butterfly" vertices that do not lie on the edge.
			EdgeMap edges = new EdgeMap ();
			int[] indices = new int[3];
			for (int i = 0; i < numFaces; ++i) {
				indices[0] = faces[i * 3];
				indices[1] = faces[i * 3 + 1];
				indices[2] = faces[i * 3 + 2];
				for (int k = 0; k < 3; ++k) {
					int index0 = Mathf.Min (indices[k], indices[(k + 1) % 3]);
					int index1 = Mathf.Max (indices[k], indices[(k + 1) % 3]);
					Vector2i edge;
					edge.x = index0;
					edge.y = index1;
					Butterfly butt = null;
					if (!edges.TryGetValue (edge, out butt)) {
						butt = new Butterfly (-1, -1, 0);
						edges.Add (edge, butt);
					}

					if (index0 == indices[k]) {
						butt.wingverts.x = indices[(k + 2) % 3]; // butterfly vert on left side of edge, ccw
						butt.count++;
					} else {
						butt.wingverts.y = indices[(k + 2) % 3]; // and right side 
						butt.count++;
					}
				}
			}

			int skipped = 0;

			TripletMap tripletMap = new TripletMap ();
			int edgeIndex = 0;
			int[] vertIdx = new int[4];
			Vector3d[] butterflyVerts = new Vector3d[4];
			//Debug.Log ("edges = " + edges.Count);
			foreach (var it in edges) {

				if (it.Value.count != 2) {
					++edgeIndex;
					continue; // not an interior edge, ignore
				}

				if (it.Value.wingverts.x < 0 || it.Value.wingverts.y < 0) {
					++edgeIndex;
					continue; // duplicate face, ignore
				}

				vertIdx[0] = it.Key.x;
				vertIdx[1] = it.Key.y;
				vertIdx[2] = it.Value.wingverts.x;
				vertIdx[3] = it.Value.wingverts.y;

				for (int i = 0; i < 4; ++i) {
					Vector3 v = verts[vertIdx[i]];
					butterflyVerts[i] = new Vector3d (v.x, v.y, v.z);
				}

				/*Matrix24*/
				MathNet.Numerics.LinearAlgebra.Matrix<double> GD = MathNet.Numerics.LinearAlgebra.Double.Matrix.Build.Dense (2, 4);
				if (!ButterflyGradDiff (butterflyVerts, ref GD)) {
					//Debug.Log(GD);
					skipped++;
					++edgeIndex;
					continue;
				}

				/*Matrix44*/
				MathNet.Numerics.LinearAlgebra.Matrix<double> GDtGD = GD.Transpose () * GD; // units will now be [m^2]
				//Debug.Log(GDtGD);

				// scatter GDtGD:
				Vector2i key = new Vector2i ();
				for (int i = 0; i < 4; i++) {
					for (int j = 0; j < 4; j++) {
						key.x = vertIdx[i];
						key.y = vertIdx[j];
						if (!tripletMap.ContainsKey (key))
							tripletMap.Add (key, 0.0);
						tripletMap[key] += GDtGD[i, j];
					}
				}
				++edgeIndex;
			}

			System.Tuple<int, int, double>[] triplets = new System.Tuple<int, int, double>[tripletMap.Count];
			int num = 0;
			foreach (var it in tripletMap) {
				triplets[num] = new System.Tuple<int, int, double> (it.Key.x, it.Key.y, it.Value);
				num++;
			}
			regularizationMatrix = MathNet.Numerics.LinearAlgebra.Double.SparseMatrix.OfIndexed ((int) numVerts, (int) numVerts, triplets);

#if DEBUG
			if (skipped > 0) {
				Debug.Log ("EdgeBasedRegularizer: skipped " + skipped + " edges out of " + edges.Count);
			}
#endif
		}

		static void DoFilterLeastSquares (
			Scene scene,
			int[] numSamplesPerInstance,
			Samples samples,
			float[] values,
			float regularizationWeight,
			float[][] vertexSignals,
			bool parallel = false
		) {
			MathNet.Numerics.Control.UseMultiThreading ();
			double massMatrixTime = 0.0, regularizationMatrixTime = 0.0, decomposeTime = 0.0, solveTime = 0.0;

			int[] sampleOffsetPerInstance = new int[scene.numInstances]; {
				int sampleOffset = 0;
				for (int i = 0; i < scene.numInstances; ++i) {
					sampleOffsetPerInstance[i] = sampleOffset;
					sampleOffset += numSamplesPerInstance[i];
				}
			}

			//#pragma omp parallel for
			if (parallel) {
				System.Threading.Tasks.Parallel.For (0, scene.numMeshes, (meshIdx) => {

					// Build reg. matrix once, it does not depend on rigid transformation matrix per instance
					MathNet.Numerics.LinearAlgebra.Double.SparseMatrix regularizationMatrix = null;
					if (regularizationWeight > 0.0f) {
						Timer regularizationMatrixTimer = new Timer ();
						BuildRegularizationMatrix (scene.meshes[meshIdx], ref regularizationMatrix, regularizationMatrixTimer);
						regularizationMatrixTime += regularizationMatrixTimer.Elapsed.TotalMilliseconds;
					}

					// Filter all the instances that point to this mesh
					//#pragma omp parallel for
					for (int i = 0; i < (int) (scene.numInstances); ++i) {
						if (scene.instances[i].meshIndex == meshIdx) {
							Timer massMatrixTimer = new Timer ();
							Timer decomposeTimer = new Timer ();
							Timer solveTimer = new Timer ();
							int sampleOffset = sampleOffsetPerInstance[i];

							// Point to samples for this instance
							Samples instanceSamples = new Samples ();
							instanceSamples.numSamples = numSamplesPerInstance[i];
							instanceSamples.samplePositions = samples.samplePositions;
							instanceSamples.sampleNormals = samples.sampleNormals;
							instanceSamples.sampleFaceNormals = samples.sampleFaceNormals;
							instanceSamples.sampleInfos = samples.sampleInfos;
							instanceSamples.sampleOffset = sampleOffset;

							FilterMeshLeastSquares (scene.meshes[meshIdx], instanceSamples, values, regularizationWeight, regularizationMatrix,
								vertexSignals[i], massMatrixTimer, decomposeTimer, solveTimer, parallel);

							massMatrixTime += massMatrixTimer.Elapsed.TotalMilliseconds;
							decomposeTime += decomposeTimer.Elapsed.TotalMilliseconds;
							solveTime += solveTimer.Elapsed.TotalMilliseconds;
						}
					}
				});

			} else {
				Timer massMatrixTimer = new Timer ();
				Timer decomposeTimer = new Timer ();
				Timer solveTimer = new Timer ();
				for (int meshIdx = 0; meshIdx < scene.numMeshes; meshIdx++) {

					// Build reg. matrix once, it does not depend on rigid transformation matrix per instance
					MathNet.Numerics.LinearAlgebra.Double.SparseMatrix regularizationMatrix = null;
					if (regularizationWeight > 0.0f) {
						Timer regularizationMatrixTimer = new Timer ();
						BuildRegularizationMatrix (scene.meshes[meshIdx], ref regularizationMatrix, regularizationMatrixTimer);
						regularizationMatrixTime = regularizationMatrixTimer.Elapsed.TotalMilliseconds;
					}

					// Filter all the instances that point to this mesh
					//#pragma omp parallel for
					for (int i = 0; i < (int) (scene.numInstances); ++i) {
						if (scene.instances[i].meshIndex == meshIdx) {
							int sampleOffset = sampleOffsetPerInstance[i];

							// Point to samples for this instance
							Samples instanceSamples = new Samples ();
							instanceSamples.numSamples = numSamplesPerInstance[i];
							instanceSamples.samplePositions = samples.samplePositions;
							instanceSamples.sampleNormals = samples.sampleNormals;
							instanceSamples.sampleFaceNormals = samples.sampleFaceNormals;
							instanceSamples.sampleInfos = samples.sampleInfos;
							instanceSamples.sampleOffset = sampleOffset;

							FilterMeshLeastSquares (scene.meshes[meshIdx], instanceSamples, values, regularizationWeight, regularizationMatrix,
								vertexSignals[i], massMatrixTimer, decomposeTimer, solveTimer);
						}
					}
				}
				massMatrixTime = massMatrixTimer.Elapsed.TotalMilliseconds;
				decomposeTime = decomposeTimer.Elapsed.TotalMilliseconds;
				solveTime = solveTimer.Elapsed.TotalMilliseconds;
			}
#if DEBUG
			Debug.Log ("Build mass matrices ...           " + massMatrixTime);
			if (regularizationWeight > 0.0f) {
				Debug.Log ("Build regularization matrices ... " + regularizationMatrixTime);
			}
			Debug.Log ("Decompose matrices ...            " + decomposeTime);
			Debug.Log ("Solve linear systems ...         " + solveTime);
#endif
		}

		static void FilterMeshLeastSquares (
			Mesh mesh,
			Samples samples,
			float[] values,
			float regularizationWeight,
			MathNet.Numerics.LinearAlgebra.Double.SparseMatrix regularization_matrix,
			float[] vertexSignals,
			Timer massMatrixTimer,
			Timer decomposeTimer,
			Timer solveTimer,
			bool parallel = false
		) {
			for (int i = 0; i < mesh.numVertices; i++)
				vertexSignals[i] = 0.0f;

#if DEBUG && UNITY_EDITOR
			if (!parallel)
				UnityEditor.EditorUtility.DisplayProgressBar ("BelfegnarInc.Baking", "Building mass_matrix ...", 0);
#endif
			massMatrixTimer.Start ();

			int[] triangles = mesh.triangles;
			TripletMap tripletMap = new TripletMap ();
			var ones = MathNet.Numerics.LinearAlgebra.Double.Vector.Build.Dense (mesh.numVertices, 1.0);

			Vector3i tri = new Vector3i ();
			for (int i = 0; i < samples.numSamples; ++i) {
				SampleInfo info = samples.sampleInfos[i + samples.sampleOffset];
				tri.x = triangles[info.triIdx * 3];
				tri.y = triangles[info.triIdx * 3 + 1];
				tri.z = triangles[info.triIdx * 3 + 2];

				float val = values[i + samples.sampleOffset] * info.dA;

				vertexSignals[tri.x] += info.bary.x * val;
				vertexSignals[tri.y] += info.bary.y * val;
				vertexSignals[tri.z] += info.bary.z * val;

				// Note: the reference paper suggests computing the mass matrix analytically.
				// Building it from samples gave smoother results for low numbers of samples per face.
				Vector2i pair;
				pair.x = tri.x;
				pair.y = tri.x;
				if (!tripletMap.ContainsKey (pair)) {
					tripletMap.Add (pair, 0.0);
				}
				tripletMap[pair] += (double) (info.bary.x * info.bary.x * info.dA);

				pair.x = tri.y;
				pair.y = tri.y;
				if (!tripletMap.ContainsKey (pair)) {
					tripletMap.Add (pair, 0.0);
				}
				tripletMap[pair] += (double) (info.bary.y * info.bary.y * info.dA);

				pair.x = tri.z;
				pair.y = tri.z;
				if (!tripletMap.ContainsKey (pair)) {
					tripletMap.Add (pair, 0.0);
				}
				tripletMap[pair] += (double) (info.bary.z * info.bary.z * info.dA);

				{
					double elem = (double) (info.bary.x * info.bary.y * info.dA);

					pair.x = tri.x;
					pair.y = tri.y;
					if (!tripletMap.ContainsKey (pair)) {
						tripletMap.Add (pair, 0.0);
					}
					tripletMap[pair] += elem;

					pair.x = tri.y;
					pair.y = tri.x;
					if (!tripletMap.ContainsKey (pair)) {
						tripletMap.Add (pair, 0.0);
					}
					tripletMap[pair] += elem;
				}

				{
					double elem = (double) (info.bary.y * info.bary.z * info.dA);

					pair.x = tri.y;
					pair.y = tri.z;
					if (!tripletMap.ContainsKey (pair)) {
						tripletMap.Add (pair, 0.0);
					}
					tripletMap[pair] += elem;

					pair.x = tri.z;
					pair.y = tri.y;
					if (!tripletMap.ContainsKey (pair)) {
						tripletMap.Add (pair, 0.0);
					}
					tripletMap[pair] += elem;
				}

				{
					double elem = (double) (info.bary.z * info.bary.x * info.dA);

					pair.x = tri.x;
					pair.y = tri.z;
					if (!tripletMap.ContainsKey (pair)) {
						tripletMap.Add (pair, 0.0);
					}
					tripletMap[pair] += elem;

					pair.x = tri.z;
					pair.y = tri.x;
					if (!tripletMap.ContainsKey (pair)) {
						tripletMap.Add (pair, 0.0);
					}
					tripletMap[pair] += elem;
				}

			}

			System.Tuple<int, int, double>[] triplets = new System.Tuple<int, int, double>[tripletMap.Count];
			int num = 0;
			foreach (var it in tripletMap) {
				triplets[num] = new System.Tuple<int, int, double> (it.Key.x, it.Key.y, it.Value);
				num++;
			}

			// Mass matrix

			//DenseMatrix mass_matrix = MathNet.Numerics.LinearAlgebra.Double.DenseMatrix.OfIndexed ((int)mesh.num_vertices, (int)mesh.num_vertices, triplets);
			MathNet.Numerics.LinearAlgebra.Double.SparseMatrix massMatrix = MathNet.Numerics.LinearAlgebra.Double.SparseMatrix.OfIndexed ((int) mesh.numVertices, (int) mesh.numVertices, triplets);

			// Fix missing data due to unreferenced verts
			{
				var lumped = massMatrix * ones;
				for (int i = 0; i < mesh.numVertices; ++i) {
					if (lumped[i] <= 0.0) { // all valid entries in mass matrix are > 0
						massMatrix[i, i] = 1.0;
					}
				}
			}

			massMatrixTimer.Stop ();

			//var order = CSparse.ColumnOrdering.MinimumDegreeAtPlusA;
			MathNet.Numerics.LinearAlgebra.Double.Factorization.SparseLDL solver = null;
			//MathNet.Numerics.LinearAlgebra.Factorization.Cholesky<double> solver = null;
			//CustomCholesky solver = null;
			//LLT solver = null;
#if DEBUG && UNITY_EDITOR
			if (!parallel)
				UnityEditor.EditorUtility.DisplayProgressBar ("BelfegnarInc.Baking", "Decompose A ...", 0.5f);
#endif
			// Optional edge-based regularization for smoother result, see paper for details
			if (regularizationWeight > 0.0f) {
				decomposeTimer.Start ();
				massMatrix.Add (regularization_matrix * regularizationWeight, massMatrix);
				solver = massMatrix.SparseLDL ();
				decomposeTimer.Stop ();
			} else {
				decomposeTimer.Start ();
				solver = massMatrix.SparseLDL ();
				decomposeTimer.Stop ();
			}

#if DEBUG && UNITY_EDITOR
			if (!parallel)
				UnityEditor.EditorUtility.DisplayProgressBar ("BelfegnarInc.Baking", "Solve ...", 0.5f);
#endif
			solveTimer.Start ();

			//Debug.Assert( solver.info() == 0/*Eigen::Success*/ );

			var b = new double[mesh.numVertices]; // Vector.Build.Dense( mesh.num_vertices );
			var x = new double[mesh.numVertices]; //Vector.Build.Dense( mesh.num_vertices );
			for (int k = 0; k < mesh.numVertices; ++k) {
				b[k] = vertexSignals[k];
				//x[k] = 0.0;
			}

			solver.Solve (b, x);

			solveTimer.Stop ();

			for (int k = 0; k < mesh.numVertices; ++k) {
				vertexSignals[k] = (float) (x[k]); // Note: allow out-of-range values
			}
#if DEBUG && UNITY_EDITOR
			if (!parallel)
				UnityEditor.EditorUtility.ClearProgressBar ();
#endif
		}

		#endregion
	}

	class Butterfly {
		public Vector2i wingverts = new Vector2i (-1, -1);
		public int count = 0;
		public Butterfly (int x = -1, int y = -1, int c = 0) {
			wingverts.x = x;
			wingverts.y = y;
			count = c;
		}
	}

}