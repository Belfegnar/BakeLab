using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Rendering;
using Timer = System.Diagnostics.Stopwatch;

namespace BelfegnarInc.BakeLab {
	public partial class Baking {
		public enum PackStyle {
			SH,
			Dir1 = 1,
			Dir3 = 3,
			Dir4
		}

		public static void ComputeIrradiance (
			Scene scene,
			Samples samples,
			PackStyle packStyle,
			int raysPerSample,
			float sceneOffset,
			float sceneMaxdistance,
			bool cpuMode,
			IList<float[]> irradValues
		) {
			// if (packStyle == PackStyle.SH)
			// 	Irradiance.BakeIrradiance (scene, samples, raysPerSample, sceneOffset, sceneMaxdistance, cpuMode, irradValues);
			// else if (packStyle == PackStyle.Dir3)
			// 	Irradiance.BakeIrradianceThreeDir (scene, samples, raysPerSample, sceneOffset, sceneMaxdistance, cpuMode, irradValues);
			// else if (packStyle == PackStyle.Dir1)
			Irradiance.BakeIrradianceOneDir (scene, samples, raysPerSample, sceneOffset, sceneMaxdistance, cpuMode, irradValues);
		}

		public static class Irradiance {

			public static void BakeIrradiance (
				Scene scene,
				Samples samples,
				int raysPerSample,
				float sceneOffset,
				float sceneMaxdistance,
				bool cpuMode,
				IList<float[]> irradValues) {
#if DEBUG
				Timer setupTimer = new Timer ();
				Timer queryTimer = new Timer ();
				setupTimer.Start ();
#endif
				var oldProbes = LightmapSettings.lightProbes;
				var oldProbeGroups = GameObject.FindObjectsOfType<LightProbeGroup> ();
				foreach (var group in oldProbeGroups)
					group.enabled = false;

				var probeGroup = (new GameObject ("BakeLabProbes")).AddComponent<LightProbeGroup> ();

				Vector3[] sampleNormals = samples.sampleNormals;
				Vector3[] samplePositions = samples.samplePositions;
				int numSamples = samples.numSamples;
				int sampleOffset = samples.sampleOffset;

				Vector3[] probePositions = new Vector3[numSamples];
				for (int i = 0; i < numSamples; i++) {
					Vector3 sampleNorm = sampleNormals[i + sampleOffset];
					Vector3 samplePos = samplePositions[i + sampleOffset];
					probePositions[i] = samplePos + sampleNorm * sceneOffset;
				}
#if UNITY_EDITOR
				probeGroup.probePositions = probePositions;
				UnityEditor.EditorUtility.SetDirty (LightmapSettings.lightProbes);
				UnityEditor.SceneView.RepaintAll ();
				var bakeResult = UnityEditor.Lightmapping.Bake ();
#endif

#if DEBUG
				Debug.Log ("Sample_offset = " + sampleOffset + ", scene_maxdistance = " + sceneMaxdistance);
				setupTimer.Stop ();
				queryTimer.Start ();
#endif

				//var probes = LightmapSettings.lightProbes.bakedProbes;
				SphericalHarmonicsL2 sh;
				for (int i = 0; i < numSamples; i++) {
#if DEBUG && UNITY_EDITOR
					if (i % 1000 == 0 && UnityEditor.EditorUtility.DisplayCancelableProgressBar ("BelfegnarInc.Baking", string.Format ("Calc irradiance values ...    {0:D}/{1:D}", i, numSamples), i * 1f / numSamples))
						return;
#endif
					//sh = probes[i];
					LightProbes.GetInterpolatedProbe (probePositions[i], null, out sh);
					int channel = 0;
					while (channel < 3) {
						int coef = 0;
						while (coef < 4) {
							irradValues[channel * 4 + coef][i + sampleOffset] = sh[channel, coef];
							coef++;
						}
						channel++;
					}
				}
#if DEBUG
				queryTimer.Stop ();
				Debug.Log ("Setup ...           " + setupTimer.Elapsed.TotalMilliseconds + " ms");
				Debug.Log ("Accum query ...     " + queryTimer.Elapsed.TotalMilliseconds + " ms");
#endif

				GameObject.DestroyImmediate (probeGroup.gameObject);
				LightmapSettings.lightProbes = oldProbes;
				foreach (var group in oldProbeGroups)
					group.enabled = true;
				UnityEditor.EditorUtility.SetDirty (LightmapSettings.lightProbes);
				UnityEditor.SceneView.RepaintAll ();
			}

			public static void BakeIrradianceThreeDir (
				Scene scene,
				Samples samples,
				int raysPerSample,
				float sceneOffset,
				float sceneMaxdistance,
				bool cpuMode,
				IList<float[]> irradValues) {
#if DEBUG
				Timer setupTimer = new Timer ();
				Timer queryTimer = new Timer ();
				setupTimer.Start ();
#endif
				var oldProbes = LightmapSettings.lightProbes;
				var oldProbeGroups = GameObject.FindObjectsOfType<LightProbeGroup> ();
				foreach (var group in oldProbeGroups)
					group.enabled = false;

				var probeGroup = (new GameObject ("BakeLabProbes")).AddComponent<LightProbeGroup> ();

				Vector3[] sampleNormals = samples.sampleNormals;
				Vector4[] sampleTangents = samples.sampleTangents;
				Vector3[] samplePositions = samples.samplePositions;
				int numSamples = samples.numSamples;
				int sampleOffset = samples.sampleOffset;

				Vector3[] probePositions = new Vector3[numSamples + 8];
				for (int i = 0; i < numSamples; i++) {
					//Vector3 sampleNorm = sampleNormals[i + sampleOffset];
					Vector3 samplePos = samplePositions[i + sampleOffset];
					probePositions[i] = samplePos; // + sampleNorm * sceneOffset;
				}
				probePositions[numSamples] = scene.bounds.min;
				probePositions[numSamples + 1] = scene.bounds.max;
				probePositions[numSamples + 2] = scene.bounds.min + new Vector3 (scene.bounds.size.x, 0, 0);
				probePositions[numSamples + 3] = scene.bounds.min + new Vector3 (0, scene.bounds.size.y, 0);
				probePositions[numSamples + 4] = scene.bounds.min + new Vector3 (0, 0, scene.bounds.size.z);
				probePositions[numSamples + 5] = scene.bounds.min + new Vector3 (scene.bounds.size.x, scene.bounds.size.y, 0);
				probePositions[numSamples + 6] = scene.bounds.min + new Vector3 (0, scene.bounds.size.y, scene.bounds.size.z);
				probePositions[numSamples + 7] = scene.bounds.min + new Vector3 (scene.bounds.size.x, 0, scene.bounds.size.z);
#if UNITY_EDITOR
				probeGroup.probePositions = probePositions;
				UnityEditor.EditorUtility.SetDirty (LightmapSettings.lightProbes);
				UnityEditor.SceneView.RepaintAll ();
				var bakeResult = UnityEditor.Lightmapping.Bake ();
#endif

#if DEBUG
				Debug.Log ("Sample_offset = " + sampleOffset + ", scene_maxdistance = " + sceneMaxdistance);
				setupTimer.Stop ();
				queryTimer.Start ();
#endif

				//var probes = LightmapSettings.lightProbes.bakedProbes;
				Vector3[] basis = {
					new Vector3 (Mathf.Sqrt (2f / 3f), 0, 1f / Mathf.Sqrt (3f)),
					new Vector3 (-1f / Mathf.Sqrt (6f), 1f / Mathf.Sqrt (2f), 1f / Mathf.Sqrt (3f)),
					new Vector3 (-1f / Mathf.Sqrt (6f), -1f / Mathf.Sqrt (2f), 1f / Mathf.Sqrt (3f))
				};
				Vector3[] dirs = new Vector3[3];
				Color[] colors = new Color[3];
				SphericalHarmonicsL2 sh;
				Onb onb;
				for (int i = 0; i < numSamples; i++) {
#if DEBUG && UNITY_EDITOR
					if (i % 1000 == 0 && UnityEditor.EditorUtility.DisplayCancelableProgressBar ("BelfegnarInc.Baking", string.Format ("Calc irradiance values ...    {0:D}/{1:D}", i, numSamples), i * 1f / numSamples))
						return;
#endif
					Vector3 sampleNorm = sampleNormals[i + sampleOffset];
					Vector4 sampleTang = sampleTangents[i + sampleOffset];
					Vector3 samplePos = samplePositions[i + sampleOffset];
					onb = new Onb (sampleNorm, sampleTang);
					dirs[0] = basis[0];
					dirs[1] = basis[1];
					dirs[2] = basis[2];
					onb.InverseTransform (ref dirs[0]);
					onb.InverseTransform (ref dirs[1]);
					onb.InverseTransform (ref dirs[2]);
					if (i == 0) {
						Debug.Log (string.Format ("{0}, {1}, {2}", dirs[0], dirs[1], dirs[2]));
					}
					LightProbes.GetInterpolatedProbe (probePositions[i] + sampleNorm * sceneOffset, null, out sh);
					sh.Evaluate (dirs, colors);
					if (i == 0) {
						Debug.Log (string.Format ("{0}, {1}, {2}", colors[0], colors[1], colors[2]));
					}
					int color = 0;
					while (color < 3) {
						int coef = 0;
						while (coef < 3) {
							irradValues[color * 3 + coef][i + sampleOffset] = colors[color][coef];
							coef++;
						}
						color++;
					}
				}
#if DEBUG
				queryTimer.Stop ();
				Debug.Log ("Setup ...           " + setupTimer.Elapsed.TotalMilliseconds + " ms");
				Debug.Log ("Accum query ...     " + queryTimer.Elapsed.TotalMilliseconds + " ms");
#endif

				GameObject.DestroyImmediate (probeGroup.gameObject);
				LightmapSettings.lightProbes = oldProbes;
				foreach (var group in oldProbeGroups)
					group.enabled = true;
				UnityEditor.EditorUtility.SetDirty (LightmapSettings.lightProbes);
				UnityEditor.SceneView.RepaintAll ();
			}

			public static void BakeIrradianceOneDir (
				Scene scene,
				Samples samples,
				int raysPerSample,
				float sceneOffset,
				float sceneMaxdistance,
				bool cpuMode,
				IList<float[]> irradValues) {
#if DEBUG
				Timer setupTimer = new Timer ();
				Timer queryTimer = new Timer ();
				setupTimer.Start ();
#endif
				var oldProbes = LightmapSettings.lightProbes;
				var oldProbeGroups = GameObject.FindObjectsOfType<LightProbeGroup> ();
				foreach (var group in oldProbeGroups)
					group.enabled = false;

				var probeGroup = (new GameObject ("BakeLabProbes")).AddComponent<LightProbeGroup> ();

				Vector3[] sampleNormals = samples.sampleNormals;
				Vector3[] samplePositions = samples.samplePositions;
				int numSamples = samples.numSamples;
				int sampleOffset = samples.sampleOffset;

				Vector3[] probePositions = new Vector3[numSamples + 8];
				for (int i = 0; i < numSamples; i++) {
					//Vector3 sampleNorm = sampleNormals[i + sampleOffset];
					Vector3 samplePos = samplePositions[i + sampleOffset];
					probePositions[i] = samplePos; // + sampleNorm * sceneOffset;
				}
				probePositions[numSamples] = scene.bounds.min;
				probePositions[numSamples + 1] = scene.bounds.max;
				probePositions[numSamples + 2] = scene.bounds.min + new Vector3 (scene.bounds.size.x, 0, 0);
				probePositions[numSamples + 3] = scene.bounds.min + new Vector3 (0, scene.bounds.size.y, 0);
				probePositions[numSamples + 4] = scene.bounds.min + new Vector3 (0, 0, scene.bounds.size.z);
				probePositions[numSamples + 5] = scene.bounds.min + new Vector3 (scene.bounds.size.x, scene.bounds.size.y, 0);
				probePositions[numSamples + 6] = scene.bounds.min + new Vector3 (0, scene.bounds.size.y, scene.bounds.size.z);
				probePositions[numSamples + 7] = scene.bounds.min + new Vector3 (scene.bounds.size.x, 0, scene.bounds.size.z);
#if UNITY_EDITOR
				probeGroup.probePositions = probePositions;
				UnityEditor.EditorUtility.SetDirty (LightmapSettings.lightProbes);
				UnityEditor.SceneView.RepaintAll ();
				var bakeResult = UnityEditor.Lightmapping.Bake ();
#endif

#if DEBUG
				Debug.Log ("Sample_offset = " + sampleOffset + ", scene_maxdistance = " + sceneMaxdistance);
				setupTimer.Stop ();
				queryTimer.Start ();
#endif

				Vector3[] dirs = new Vector3[1];
				Color[] colors = new Color[1];
				SphericalHarmonicsL2 sh;
				for (int i = 0; i < numSamples; i++) {
#if DEBUG && UNITY_EDITOR
					if (i % 1000 == 0 && UnityEditor.EditorUtility.DisplayCancelableProgressBar ("BelfegnarInc.Baking", string.Format ("Calc irradiance values ...    {0:D}/{1:D}", i, numSamples), i * 1f / numSamples))
						return;
#endif
					Vector3 sampleNorm = sampleNormals[i + sampleOffset];
					Vector3 samplePos = samplePositions[i + sampleOffset];
					dirs[0] = sampleNorm;

					LightProbes.GetInterpolatedProbe (probePositions[i] + sampleNorm * sceneOffset, null, out sh);
					sh.Evaluate (dirs, colors);
					int coef = 0;
					while (coef < 3) {
						irradValues[coef][i + sampleOffset] = colors[0][coef];
						coef++;
					}
				}
#if DEBUG
				queryTimer.Stop ();
				Debug.Log ("Setup ...           " + setupTimer.Elapsed.TotalMilliseconds + " ms");
				Debug.Log ("Accum query ...     " + queryTimer.Elapsed.TotalMilliseconds + " ms");
#endif

				GameObject.DestroyImmediate (probeGroup.gameObject);
				LightmapSettings.lightProbes = oldProbes;
				foreach (var group in oldProbeGroups)
					group.enabled = true;
				UnityEditor.EditorUtility.SetDirty (LightmapSettings.lightProbes);
				UnityEditor.SceneView.RepaintAll ();
			}
		}
	}
}