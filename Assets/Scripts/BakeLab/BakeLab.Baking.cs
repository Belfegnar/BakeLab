using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Timer = System.Diagnostics.Stopwatch;

namespace BelfegnarInc.BakeLab {
	public static partial class Baking {

		public static void ComputeAO (
			Scene scene,
			Samples samples,
			int raysPerSample,
			float sceneOffset,
			float sceneMaxdistance,
			bool cpuMode,
			float[] aoValues
		) {
			BakeAO (scene, samples, raysPerSample, sceneOffset, sceneMaxdistance, cpuMode, aoValues);
		}

		public static void ComputeRadiance (
			Scene scene,
			Samples samples,
			int raysPerSample,
			float sceneOffset,
			float sceneMaxdistance,
			bool cpuMode,
			float[] radValues
		) {
			BakeRadiance (scene, samples, raysPerSample, sceneOffset, sceneMaxdistance, cpuMode, radValues);
		}

		public static void ComputeTangentSpaceLightDirection (
			Scene scene,
			Samples samples,
			int raysPerSample,
			float sceneOffset,
			float sceneMaxdistance,
			bool cpuMode,
			float[] dirValues,
			int attrIndex
		) {
			BakeDistToLightInTangentSpace (scene, samples, raysPerSample, sceneOffset, sceneMaxdistance, cpuMode, dirValues, attrIndex);
		}

		/// <summary>
		/// Generate random unsigned int in [0, 2^24)
		/// </summary>
		static uint Lcg (ref uint prev) {
			const uint LCG_A = 1664525u;
			const uint LCG_C = 1013904223u;
			prev = (LCG_A * prev + LCG_C);
			return prev & 0x00FFFFFF;
		}

		/// <summary>
		/// Generate random float in [0, 1)
		/// </summary>
		static float Rnd (ref uint prev) {
			return ((float) Lcg (ref prev) / (float) 0x01000000);
		}

		static void CosineSampleHemisphere (float u1, float u2, ref Vector3 p) {
			// Uniformly sample disk.
			float r = Mathf.Sqrt (u1);
			float phi = 2.0f * Mathf.PI * u2;
			p.x = r * Mathf.Cos (phi);
			p.y = r * Mathf.Sin (phi);

			// Project up to hemisphere.
			p.z = Mathf.Sqrt (Mathf.Max (0.0f, 1.0f - p.x * p.x - p.y * p.y));
		}

		static void BakeAO (
			Scene scene,
			Samples samples,
			int raysPerSample,
			float sceneOffset,
			float sceneMaxdistance,
			bool cpuMode,
			float[] aoValues) {
#if DEBUG
			Timer setupTimer = new Timer ();
			Timer queryTimer = new Timer ();
			setupTimer.Start ();
#endif

			uint seed = 0;
			int sqrtRaysPerSample = (int) (Mathf.Sqrt ((float) (raysPerSample)) + .5f);
			Vector3[] rotations = new Vector3[sqrtRaysPerSample * sqrtRaysPerSample];
			int idx = 0;
			Onb onb = new Onb (Vector3.forward);
			//			Transform onbGO = (new GameObject ("Onb")).transform;

			for (int i = 0; i < sqrtRaysPerSample; ++i)
				for (int j = 0; j < sqrtRaysPerSample; ++j) {
					Vector3 rayDir;
					rayDir.x = 0;
					rayDir.y = 0;
					rayDir.z = 0;
					float u0 = ((float) (i) + Rnd (ref seed)) / (float) (sqrtRaysPerSample);
					float u1 = ((float) (j) + Rnd (ref seed)) / (float) (sqrtRaysPerSample);

					CosineSampleHemisphere (u0, u1, ref rayDir);
					u0 = Rnd (ref seed);
					u1 = Rnd (ref seed);
					rotations[idx++] = rayDir; //Random.onUnitSphere;
				}

			Vector3[] sampleNormals = samples.sampleNormals;
			Vector4[] sampleTangents = samples.sampleTangents;
			Vector3[] samplePositions = samples.samplePositions;
			int numSamples = samples.numSamples;
			int sampleOffset = samples.sampleOffset;

#if DEBUG
			Debug.Log ("Sample_offset = " + sampleOffset + ", scene_maxdistance = " + sceneMaxdistance);
			setupTimer.Stop ();
			queryTimer.Start ();
#endif

			Vector3 dir = Vector3.zero;
			UnityEngine.Ray ray = new UnityEngine.Ray ();
			int wrongDirections = 0;
			for (int i = 0; i < numSamples; i++) {
#if DEBUG && UNITY_EDITOR
				if (i % 1000 == 0 && UnityEditor.EditorUtility.DisplayCancelableProgressBar ("BelfegnarInc.Baking", string.Format ("Calc rays ...    {0:D}/{1:D}", i, numSamples), i * 1f / numSamples))
					return;
#endif
				float hits = 0;
				for (int s = 0; s < raysPerSample; s++) {
					dir = rotations[s];
					Vector3 sampleNorm = sampleNormals[i + sampleOffset];
					Vector4 sampleTang = sampleTangents[i + sampleOffset];
					Vector3 samplePos = samplePositions[i + sampleOffset];

					onb = new Onb (sampleNorm, sampleTang);
					onb.InverseTransform (ref dir);

					//					onbGO.forward = sampleNorm;
					//					dir = onbGO.TransformDirection (dir);

#if false
					ray.origin = samplePos + sampleNorm * sceneOffset + sceneMaxdistance * dir;
					ray.direction = -dir;

					if (Physics.Raycast (ray, sceneMaxdistance - sceneOffset)) {
						hits += 1f;
					}
#else
					ray.origin = samplePos + sampleNorm * sceneOffset;
					ray.direction = dir;

					if (Vector3.Dot (dir, sampleNorm) < 0) {
						wrongDirections++;
						//Debug.DrawLine (Vector3.zero, dir, Color.red, 90, false);
					} else
					if (Physics.Raycast (ray, sceneMaxdistance)) {
						hits += 1f;
					}

#endif
				}
				aoValues[i + sampleOffset] = 1.0f - hits / raysPerSample;
			}
#if DEBUG
			queryTimer.Stop ();
			Debug.Log ("WRONG_DIRECTIONS: " + wrongDirections);
			Debug.Log ("Setup ...           " + setupTimer.Elapsed.TotalMilliseconds + " ms");
			Debug.Log ("Accum query ...     " + queryTimer.Elapsed.TotalMilliseconds + " ms");
#endif
		}

		static void BakeRadiance (
			Scene scene,
			Samples samples,
			int raysPerSample,
			float sceneOffset,
			float sceneMaxdistance,
			bool cpuMode,
			float[] radValues) {
#if DEBUG
			Timer setupTimer = new Timer ();
			Timer queryTimer = new Timer ();
			setupTimer.Start ();
#endif

			var lights = GameObject.FindObjectsOfType<Light> ();
			Vector3 lightDir = Vector3.zero;
			float intensity = 0;
			for (int l = 0; l < lights.Length; l++) {
				if (lights[l].type == LightType.Directional) {
					intensity = lights[l].intensity * lights[l].color.grayscale;
					lightDir = lights[l].transform.forward;
					for (int l2 = l + 1; l2 < lights.Length; l2++) {
						if (lights[l2].type == LightType.Directional) {
							float intensity2 = lights[l2].intensity * lights[l2].color.grayscale;
							if (intensity2 > intensity) {
								intensity = intensity2;
								lightDir = lights[l2].transform.forward;
							}
						}
					}
					break;
				}
			}

			Vector3[] sampleNormals = samples.sampleNormals;
			Vector3[] sampleFaceNormals = samples.sampleFaceNormals;
			Vector3[] samplePositions = samples.samplePositions;
			int numSamples = samples.numSamples;
			int sampleOffset = samples.sampleOffset;
			Vector3 dir = -lightDir;
			UnityEngine.Ray ray = new UnityEngine.Ray ();

#if DEBUG
			Debug.Log ("Sample_offset = " + sampleOffset + ", scene_maxdistance = " + sceneMaxdistance);
			setupTimer.Stop ();
			queryTimer.Start ();
#endif

			for (int i = 0; i < numSamples; i++) {
#if DEBUG && UNITY_EDITOR
				if (i % 1000 == 0 && UnityEditor.EditorUtility.DisplayCancelableProgressBar ("BelfegnarInc.Baking", string.Format ("Calc rays ...    {0:D}/{1:D}", i, numSamples), i * 1f / numSamples))
					return;
#endif

				Vector3 sampleNorm = sampleNormals[i + sampleOffset];
				Vector3 sampleFaceNorm = sampleFaceNormals[i + sampleOffset];
				Vector3 samplePos = samplePositions[i + sampleOffset];
				float hits = 0.0f;
				ray.origin = samplePos + sampleFaceNorm * sceneOffset; // + dir*scene_maxdistance;
				ray.direction = dir;
				if (Vector3.Dot (sampleNorm, dir) < 0 || Physics.Raycast (ray, sceneMaxdistance)) {
					hits += 1f;
				}
				radValues[i + sampleOffset] = (1.0f - hits) * intensity;
			}
#if DEBUG
			queryTimer.Stop ();
			Debug.Log ("Setup ...           " + setupTimer.Elapsed.TotalMilliseconds + " ms");
			Debug.Log ("Accum query ...     " + queryTimer.Elapsed.TotalMilliseconds + " ms");
#endif
		}

		static void BakeDistToLightInTangentSpace (
			Scene scene,
			Samples samples,
			int raysPerSample,
			float sceneOffset,
			float sceneMaxdistance,
			bool cpuMode,
			float[] dirValues,
			int attrIndex) {
#if DEBUG
			Timer setupTimer = new Timer ();
			Timer queryTimer = new Timer ();
			setupTimer.Start ();
#endif

			int numPointLights = 0;
			var lights = GameObject.FindObjectsOfType<Light> ();
			Vector3 lightDir = Vector3.zero;
			float intensity = 0;
			for (int l = 0; l < lights.Length; l++) {
				if (lights[l].type == LightType.Directional) {
					intensity = lights[l].intensity; // * lights[l].color.grayscale;
					lightDir = lights[l].transform.forward;
					for (int l2 = l + 1; l2 < lights.Length; l2++) {
						if (lights[l2].type == LightType.Directional) {
							float intensity2 = lights[l2].intensity; // * lights[l2].color.grayscale;
							if (intensity2 > intensity) {
								intensity = intensity2;
								lightDir = lights[l2].transform.forward;
							}
						}
					}
					break;
				}
			}

			for (int l = 0; l < lights.Length; l++) {
				if (lights[l].type == LightType.Point) {
					numPointLights++;
				}
			}
			Vector4[] lightsPositions = new Vector4[numPointLights];
			int cL = 0;
			for (int l = 0; l < lights.Length; l++) {
				if (lights[l].type == LightType.Point) {
					lightsPositions[l] = lights[l].transform.position;
					lightsPositions[l].w = lights[l].range;
					cL++;
				}
			}

			//			Transform onbGO = (new GameObject ("Onb")).transform;
			Vector3[] sampleNormals = samples.sampleNormals;
			Vector4[] sampleTangents = samples.sampleTangents;
			Vector3[] samplePositions = samples.samplePositions;
			int numSamples = samples.numSamples;
			int sampleOffset = samples.sampleOffset;
			Vector3 dir = lightDir;

#if DEBUG
			Debug.Log ("dir = " + (lightDir * intensity) + ", intensity = " + intensity);
			Debug.Log ("LightDir = " + lightDir + ", sample_offset = " + sampleOffset + ", scene_maxdistance = " + sceneMaxdistance);
			setupTimer.Stop ();
			queryTimer.Start ();
#endif

			for (int i = 0; i < numSamples; i++) {
#if DEBUG && UNITY_EDITOR
				if (i % 1000 == 0 && UnityEditor.EditorUtility.DisplayCancelableProgressBar ("BelfegnarInc.Baking", string.Format ("Calc rays ...    {0:D}/{1:D}", i, numSamples), i * 1f / numSamples))
					return;
#endif

				Vector3 sampleNorm = sampleNormals[i + sampleOffset];
				Vector4 sampleTang = sampleTangents[i + sampleOffset];
				Vector3 samplePos = samplePositions[i + sampleOffset];

				dir = lightDir * intensity;
				Vector4 lightPos = samplePos - lightDir;
				float minDist = float.MaxValue;
				float inten = 0;
				for (int l = 0; l < numPointLights; l++) {
					Vector4 lp = lightsPositions[l];
					float dist = (samplePos.x - lp.x) * (samplePos.x - lp.x) + (samplePos.y - lp.y) * (samplePos.y - lp.y) + (samplePos.z - lp.z) * (samplePos.z - lp.z);
					if (dist < lp.w * lp.w && dist < minDist) {
						minDist = dist;
						lightPos = lp;
						inten = 1f - Mathf.Sqrt (dist) / lp.w;
					}
				}

				if (inten <= 0 && numPointLights != 0) {
					dirValues[i + sampleOffset] = 0;
				} else {
					float mag = 1f;
					if (numPointLights != 0) {
						dir += (samplePos - lightPos.xyz ()).normalized * inten;
						mag = dir.magnitude;
						dir *= 1f / mag; //.Normalize ();
					}

					Onb onb = new Onb (sampleNorm, sampleTang);
					onb.Transform (ref dir);

					dirValues[i + sampleOffset] = dir[attrIndex] * mag;
				}
			}
#if DEBUG
			queryTimer.Stop ();
			Debug.Log ("Setup ...           " + setupTimer.Elapsed.TotalMilliseconds + " ms");
			Debug.Log ("Accum query ...     " + queryTimer.Elapsed.TotalMilliseconds + " ms");
#endif
		}

		public static Vector3 xyz (this Vector4 v) {
			Vector3 result;
			result.x = v.x;
			result.y = v.y;
			result.z = v.z;
			return result;
		}
	}

	/// <summary>
	/// Orthonormal basis
	/// </summary>
	struct Onb {

		Vector3 _tangent;
		Vector3 _binormal;
		Vector3 _normal;

		public Onb (Vector3 normal) {
			_normal = normal;

			if (Mathf.Abs (_normal.x) > Mathf.Abs (_normal.z)) {
				_binormal.x = -_normal.y;
				_binormal.y = _normal.x;
				_binormal.z = 0;
			} else {
				_binormal.x = 0;
				_binormal.y = -_normal.z;
				_binormal.z = _normal.y;
			}

			_binormal.Normalize ();
			_tangent = Vector3.Cross (_binormal, _normal);
		}

		public Onb (Vector3 normal, Vector4 tangent) {
			_normal = normal;
			_tangent.x = tangent.x;
			_tangent.y = tangent.y;
			_tangent.z = tangent.z;
			_binormal = Vector3.Cross (normal, _tangent) * tangent.w;
		}

		public Onb (Vector3 normal, Vector3 tangent, float tangentW) {
			_normal = normal;
			_tangent.x = tangent.x;
			_tangent.y = tangent.y;
			_tangent.z = tangent.z;
			_binormal = Vector3.Cross (normal, _tangent) * tangentW;
		}

		public void Transform (ref Vector3 p) {
			Vector3 result;
			result.x = p.x * _tangent.x + p.y * _tangent.y + p.z * _tangent.z;
			result.y = p.x * _binormal.x + p.y * _binormal.y + p.z * _binormal.z;
			result.z = p.x * _normal.x + p.y * _normal.y + p.z * _normal.z;
			p = result;
		}

		public void InverseTransform (ref Vector3 p) {
			Vector3 result;
			result.x = p.x * _tangent.x + p.y * _binormal.x + p.z * _normal.x;
			result.y = p.x * _tangent.y + p.y * _binormal.y + p.z * _normal.y;
			result.z = p.x * _tangent.z + p.y * _binormal.z + p.z * _normal.z;
			p = result;
			//p = p.x*_tangent + p.y*_binormal + p.z*_normal;
		}
	}

}