using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace BelfegnarInc.BakeLab {
	/// <summary>
	/// Contains scene sampling methods.
	/// </summary>
	public static class Sample {

		/// <summary>
		/// Dispributes samples on scene
		/// </summary>
		/// <param name="scene">Scene object</param>
		/// <param name="minSamplesPerTriangle">Number of samples per-triangle needed</param>
		/// <param name="requestedNumSamples">Number of samples needed in total</param>
		/// <param name="numSamplesPerInstance">Output array of numbers of samples per-instance</param>
		public static int DistributeSamples (
			Scene scene,
			int minSamplesPerTriangle,
			int requestedNumSamples,
			int[] numSamplesPerInstance
		) {
			return Distribute (scene, minSamplesPerTriangle, requestedNumSamples, numSamplesPerInstance);
		}

		/// <summary>
		/// Generates sample points for instance
		/// </summary>
		/// <param name="scene">Scene object</param>
		/// <param name="numSamplesPerInstance">Array of numbers of samples per-instance</param>
		/// <param name="minSamplesPerTriangle">Number of samples per-triangle needed</param>
		/// <param name="samples">Output Samples object</param>
		public static void SampleInstances (
			Scene scene,
			int[] numSamplesPerInstance,
			int minSamplesPerTriangle,
			Samples samples
		) {
			SampleInstances (scene, minSamplesPerTriangle, samples, numSamplesPerInstance);
		}

		// Ref: https://en.wikipedia.org/wiki/Halton_sequence
		static float Halton2 (uint index) {
			float result = 0.0f;
			const float invBase = 1.0f / 2;
			float f = invBase;
			uint i = index;
			while (i > 0) {
				result += f * (i % 2);
				i = i / 2;
				f *= invBase;
			}
			return result;
		}

		static float Halton3 (uint index) {
			float result = 0.0f;
			const float invBase = 1.0f / 3;
			float f = invBase;
			uint i = index;
			while (i > 0) {
				result += f * (i % 3);
				i = i / 3;
				f *= invBase;
			}
			return result;
		}

		static uint Tea2 (uint val0, uint val1) {
			uint v0 = val0;
			uint v1 = val1;
			uint s0 = 0;

			for (uint n = 0; n < 2; n++) {
				s0 += 0x9e3779b9;
				v0 += ((v1 << 4) + 0xa341316c) ^ (v1 + s0) ^ ((v1 >> 5) + 0xc8013ea4);
				v1 += ((v0 << 4) + 0xad90777d) ^ (v0 + s0) ^ ((v0 >> 5) + 0x7e95761e);
			}

			return v0;
		}

		static uint Tea4 (uint val0, uint val1) {
			uint v0 = val0;
			uint v1 = val1;
			uint s0 = 0;

			for (uint n = 0; n < 4; n++) {
				s0 += 0x9e3779b9;
				v0 += ((v1 << 4) + 0xa341316c) ^ (v1 + s0) ^ ((v1 >> 5) + 0xc8013ea4);
				v1 += ((v0 << 4) + 0xad90777d) ^ (v0 + s0) ^ ((v0 >> 5) + 0x7e95761e);
			}

			return v0;
		}

		// Generate random unsigned int in [0, 2^24)
		static uint Lcg (ref uint prev) {
			const uint LCG_A = 1664525u;
			const uint LCG_C = 1013904223u;
			prev = (LCG_A * prev + LCG_C);
			return prev & 0x00FFFFFF;
		}

		// Generate random float in [0, 1)
		static float Rnd (ref uint prev) {
			return ((float) Lcg (ref prev) / (float) 0x01000000);
		}

		static Vector3 Faceforward (Vector3 normal, Vector3 geom_normal) {
			if (Vector3.Dot (normal, geom_normal) > 0.0f) return normal;
			return -normal;
		}

		static double TriangleArea (Vector3 v0, Vector3 v1, Vector3 v2) {
			Vector3 e0 = v1 - v0;
			Vector3 e1 = v2 - v0;
			Vector3 c = Vector3.Cross (e0, e1);
			double x = c.x, y = c.y, z = c.z;
			return 0.5 * System.Math.Sqrt (x * x + y * y + z * z);
		}

		static int Distribute (
			Scene scene,
			int minSamplesPerTriangle,
			int requestedNumSamples,
			int[] numSamplesPerInstance
		) {
			// Compute min samples per instance
			int[] minSamplesPerInstance = new int[scene.numInstances];
			int numTriangles = 0;
			for (int i = 0; i < scene.numInstances; ++i) {
				Mesh mesh = scene.meshes[scene.instances[i].meshIndex];
				minSamplesPerInstance[i] = (minSamplesPerTriangle * mesh.numTriangles);
				numTriangles += mesh.numTriangles;
			}
			int minNumSamples = minSamplesPerTriangle * numTriangles;
			int numSamples = Mathf.Max (minNumSamples, requestedNumSamples);

			// Compute surface area per instance.
			// Note: for many xforms, we could compute surface area per mesh instead of per instance.
			double[] areaPerInstance = new double[scene.numInstances];
			Vector3Int tri = new Vector3Int ();
			if (numSamples > minNumSamples) {

				//#pragma omp parallel for
				//System.Threading.Tasks.Parallel.For (0, scene.numInstances, (idx, state) =>
				for (int idx = 0; idx < scene.numInstances; ++idx) {
					Mesh mesh = scene.meshes[scene.instances[idx].meshIndex];
					Matrix4x4 matrix = scene.instances[idx].matrix;
					int[] triangles = mesh.triangles;

					for (int tri_idx = 0; tri_idx < mesh.numTriangles; ++tri_idx) {
						tri.x = triangles[tri_idx * 3];
						tri.y = triangles[tri_idx * 3 + 1];
						tri.z = triangles[tri_idx * 3 + 2];
						double area = TriangleArea (matrix.MultiplyPoint (mesh.vertices[tri.x]), matrix.MultiplyPoint (mesh.vertices[tri.y]), matrix.MultiplyPoint (mesh.vertices[tri.z]));
						areaPerInstance[idx] += area;
					}
				} //);
			}

			// Distribute samples
			InstanceSamplerCallback cb = new InstanceSamplerCallback (minSamplesPerInstance, areaPerInstance);
			DistributeGeneric (cb, numSamples, scene.numInstances, numSamplesPerInstance);

			return numSamples;
		}

		/// <summary>
		/// Distributes samples per element (instance, triangle, ...), satisfying two constraints:
		///  1) total samples add up to numSamples.
		///  2) minimum number of samples per element.
		/// Any extra samples are placed according to area ratios.
		/// The template param is for specifying area per element, and min samples per element.
		/// </summary>
		/// <param name="callback">Distribution callback</param>
		/// <param name="numSamples">Number of samples</param>
		/// <param name="numElements">Number of elements</param>
		/// <param name="numSamplesPerElement">Output array of numbers of samples per-element</param>
		static void DistributeGeneric (
			SamplerCallback callback,
			int numSamples,
			int numElements,
			int[] numSamplesPerElement
		) {
			// First place minimum samples per element
			int sampleCount = 0; // For all elements
			for (int i = 0; i < numElements; ++i) {
				numSamplesPerElement[i] = callback.MinSamples (i);
				sampleCount += numSamplesPerElement[i];
			}

#if DEBUG
			Debug.Assert (numSamples >= sampleCount);
#endif

			if (numSamples > sampleCount) {

				// Area-based sampling
				int num_area_based_samples = numSamples - sampleCount;

				// Compute surface area of each element
				double total_area = 0.0;
				for (int i = 0; i < numElements; ++i) {
					total_area += callback.Area (i);
				}

				// Distribute
				for (int i = 0; i < numElements && sampleCount < numSamples; ++i) {
					int n = Mathf.Min (numSamples - sampleCount, (int) (num_area_based_samples * callback.Area (i) / total_area));
					numSamplesPerElement[i] += n;
					sampleCount += n;
				}

				// There could be a few samples left over. Place one sample per element until target sample count is reached.
#if DEBUG
				Debug.Assert (numSamples - sampleCount <= numElements);
#endif
				for (int i = 0; i < numElements && sampleCount < numSamples; ++i) {
					numSamplesPerElement[i] += 1;
					sampleCount += 1;
				}
			}

#if DEBUG
			Debug.Assert (sampleCount == numSamples);
#endif
		}

		static void SampleInstances (
			Scene scene,
			int min_samples_per_triangle,
			Samples ao_samples,
			int[] numSamplesPerInstance
		) {

			int[] sampleOffsets = new int[scene.numInstances]; {
				int sampleOffset = 0;
				for (int i = 0; i < scene.numInstances; ++i) {
					sampleOffsets[i] = sampleOffset;
					sampleOffset += numSamplesPerInstance[i];
				}
			}

			//#pragma omp parallel for
			for (int i = 0; i < scene.numInstances; ++i) {
				int sampleOffset = sampleOffsets[i];
				// Point to samples for this instance
				Samples instanceAoSamples = new Samples ();
				instanceAoSamples.numSamples = numSamplesPerInstance[i];
				instanceAoSamples.samplePositions = ao_samples.samplePositions;
				instanceAoSamples.sampleNormals = ao_samples.sampleNormals;
				instanceAoSamples.sampleTangents = ao_samples.sampleTangents;
				instanceAoSamples.sampleFaceNormals = ao_samples.sampleFaceNormals;
				instanceAoSamples.sampleInfos = ao_samples.sampleInfos;
				instanceAoSamples.sampleOffset = sampleOffset;

				Matrix4x4 matrix = scene.instances[i].matrix;
				SampleInstance (scene.meshes[scene.instances[i].meshIndex], matrix, (uint) i, min_samples_per_triangle, instanceAoSamples);
			}
		}

		static void SampleInstance (
			Mesh mesh,
			Matrix4x4 matrix,
			uint seed,
			int minSamplesPerTriangle,
			Samples samples
		) {

			// Setup access to mesh data
			Matrix4x4 matrixInvtrans = matrix.inverse.transpose;
#if DEBUG
			Debug.Assert (samples.numSamples >= mesh.numTriangles * minSamplesPerTriangle);
#endif

			int[] triangles = mesh.triangles;
			Vector3[] samplePositions = samples.samplePositions;
			Vector3[] sampleNorms = samples.sampleNormals;
			Vector4[] sampleTangs = samples.sampleTangents;
			Vector3[] sampleFaceNorms = samples.sampleFaceNormals;
			SampleInfo[] sampleInfos = samples.sampleInfos;
			int sampleOffset = samples.sampleOffset;

			// Compute triangle areas
			double[] triAreas = new double[mesh.numTriangles];
			Vector3Int tri = new Vector3Int ();
			for (int tri_idx = 0; tri_idx < mesh.numTriangles; tri_idx++) {
				tri.x = triangles[tri_idx * 3];
				tri.y = triangles[tri_idx * 3 + 1];
				tri.z = triangles[tri_idx * 3 + 2];
				double area = TriangleArea (matrix.MultiplyPoint (mesh.vertices[tri.x]), matrix.MultiplyPoint (mesh.vertices[tri.y]), matrix.MultiplyPoint (mesh.vertices[tri.z]));
				triAreas[tri_idx] = area;
			}

			// Get sample counts
			int[] triSampleCounts = new int[mesh.numTriangles];
			TriangleSamplerCallback cb = new TriangleSamplerCallback (minSamplesPerTriangle, triAreas);
			DistributeGeneric (cb, samples.numSamples, mesh.numTriangles, triSampleCounts);

			// Place samples
			int sampleIdx = 0;
			Vector3[] verts = new Vector3[3];
			Vector3[] norms = new Vector3[3];
			Vector4[] tangents = new Vector4[3];
			for (int triIdx = 0; triIdx < mesh.numTriangles; triIdx++) {
				tri.x = triangles[triIdx * 3];
				tri.y = triangles[triIdx * 3 + 1];
				tri.z = triangles[triIdx * 3 + 2];
				verts[0] = mesh.vertices[tri.x];
				verts[1] = mesh.vertices[tri.y];
				verts[2] = mesh.vertices[tri.z];
				Vector3[] normals = null;

				if (mesh.normals != null) {
					norms[0] = mesh.normals[tri.x];
					norms[1] = mesh.normals[tri.y];
					norms[2] = mesh.normals[tri.z];

					tangents[0] = mesh.tangents[tri.x];
					tangents[1] = mesh.tangents[tri.y];
					tangents[2] = mesh.tangents[tri.z];
					normals = norms;
				}
				SampleTriangle (matrix, matrixInvtrans, verts, normals, tangents,
					triIdx, triSampleCounts[triIdx], triAreas[triIdx],
					seed,
					samplePositions, sampleNorms, sampleFaceNorms, sampleTangs, sampleInfos, sampleIdx, sampleOffset);
				sampleIdx += triSampleCounts[triIdx];
			}

#if DEBUG
			Debug.Assert (sampleIdx == samples.numSamples);
#if DEBUG_MESH_SAMPLES
			for (int i = 0; i < samples.numSamples; ++i) {
				SampleInfo info = sampleInfos[i];
				Debug.Log ("sample info (" + i + "): " + info.triIdx + ", (" + info.bary[0] + ", " + info.bary[1] + ", " + info.bary[2] + "), " + info.dA);
			}
#endif
#endif
		}

		static void SampleTriangle (Matrix4x4 matrix, Matrix4x4 matrixInvtrans,
			IList<Vector3> verts, IList<Vector3> normals, IList<Vector4> tangents,
			int triIdx, int triSampleCount, double triArea,
			uint baseSeed,
			Vector3[] samplePositions, Vector3[] sampleNorms, Vector3[] sampleFaceNorms, Vector4[] sampleTangs, SampleInfo[] sampleInfos, int sampleId, int sampleOffset) {
			Vector3 v0 = verts[0];
			Vector3 v1 = verts[1];
			Vector3 v2 = verts[2];

			Vector3 faceNormal = Vector3.Normalize (Vector3.Cross (v1 - v0, v2 - v0));
			Vector3 n0, n1, n2;
			if (normals != null) {
				n0 = Faceforward (normals[0], faceNormal);
				n1 = Faceforward (normals[1], faceNormal);
				n2 = Faceforward (normals[2], faceNormal);
			} else {
				// missing vertex normals, so use face normal.
				n0 = faceNormal;
				n1 = faceNormal;
				n2 = faceNormal;
			}

			Vector4 t0 = tangents[0];
			Vector4 t1 = tangents[1];
			Vector4 t2 = tangents[2];

			// Random offset per triangle, to shift Halton points
			uint seed = Tea4 (baseSeed, (uint) triIdx);
			Vector2 offset = new Vector2 (Rnd (ref seed), Rnd (ref seed));

			for (int index = 0; index < triSampleCount; ++index) {
				sampleInfos[sampleOffset + sampleId + index].triIdx = (uint) triIdx;
				sampleInfos[sampleOffset + sampleId + index].dA = (float) (triArea / triSampleCount);

				// Random point in unit square
				float r1 = offset.x + Halton2 ((uint) index + 1);
				r1 = r1 - (int) r1;
				float r2 = offset.y + Halton3 ((uint) index + 1);
				r2 = r2 - (int) r2;
#if DEBUG
				Debug.Assert (r1 >= 0 && r1 <= 1);
				Debug.Assert (r2 >= 0 && r2 <= 1);
#endif

				// Map to triangle. Ref: PBRT 2nd edition, section 13.6.4
				Vector3 bary = sampleInfos[sampleOffset + sampleId + index].bary;
				float sqrt_r1 = Mathf.Sqrt (r1);
				bary.x = 1.0f - sqrt_r1;
				bary.y = r2 * sqrt_r1;
				bary.z = 1.0f - bary.x - bary.y;
				sampleInfos[sampleOffset + sampleId + index].bary = bary;

				samplePositions[sampleOffset + sampleId + index] = matrix.MultiplyPoint (bary.x * v0 + bary.y * v1 + bary.z * v2);
				sampleNorms[sampleOffset + sampleId + index] = Vector3.Normalize (matrixInvtrans.MultiplyVector (bary.x * n0 + bary.y * n1 + bary.z * n2));
				sampleTangs[sampleOffset + sampleId + index] = Vector3.Normalize (matrixInvtrans.MultiplyVector (bary.x * t0 + bary.y * t1 + bary.z * t2));
				sampleTangs[sampleOffset + sampleId + index].w = bary.x * t0.w + bary.y * t1.w + bary.z * t2.w;
				sampleFaceNorms[sampleOffset + sampleId + index] = Vector3.Normalize (matrixInvtrans.MultiplyVector (faceNormal));
			}
		}
	}

	public interface SamplerCallback {
		int MinSamples (int i);
		double Area (int i);
	}

	class InstanceSamplerCallback : SamplerCallback {
		public InstanceSamplerCallback (int[] minSamplesPerInstance, double[] areaPerInstance) {
			_minSamplesPerInstance = minSamplesPerInstance;
			_areaPerInstance = areaPerInstance;
		}

		public int MinSamples (int i) {
			return _minSamplesPerInstance[i];
		}
		public double Area (int i) {
			return _areaPerInstance[i];
		}

		int[] _minSamplesPerInstance;
		double[] _areaPerInstance;
	}

	class TriangleSamplerCallback : SamplerCallback {
		public TriangleSamplerCallback (int minSamplesPerTriangle, double[] areaPerTriangle) {
			_minSamplesPerTriangle = minSamplesPerTriangle;
			_areaPerTriangle = areaPerTriangle;
		}

		public int MinSamples (int i) {
			return _minSamplesPerTriangle; // same for every triangle
		}
		public double Area (int i) {
			return _areaPerTriangle[i];
		}

		int _minSamplesPerTriangle;
		double[] _areaPerTriangle;
	}

	public class Samples {
		public int numSamples;
		public Vector3[] samplePositions;
		public Vector3[] sampleNormals;
		public Vector3[] sampleFaceNormals;
		public Vector4[] sampleTangents;
		public SampleInfo[] sampleInfos;
		public int sampleOffset;
	}

	public struct SampleInfo {
		public uint triIdx;
		public Vector3 bary;
		public float dA;
	}
}