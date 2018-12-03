using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Timer = System.Diagnostics.Stopwatch;

namespace BelfegnarInc.BakeLab {
	/// <summary>
	/// Contains main methods.
	/// </summary>
	public class BakeLabMain : MonoBehaviour {
		/// <summary>
		/// Actual configuration settings.
		/// </summary>
		public BakeLabConfig config;
		public bool parallel;
		public Transform lightTransform;
		public Vector3 fromRotation;
		public Vector3 toRotation;

		/// <summary>
		/// Swap two params (Generic)
		/// </summary>
		/// <param name="a">First parameter</param>
		/// <param name="b">Second parameter</param>
		void Swap<T> (ref T a, ref T b) {
			T c = a;
			a = b;
			b = c;
		}

		/// <summary>
		/// Finds MeshFilters in child GOs, creates scene/instance data needed
		/// </summary>
		/// <param name="scene">Empty scene object</param>
		/// <param name="sceneBounds">Bounds, encaps. founded meshes</param>
		/// <returns>True, if no errors found</returns>
		bool LoadScene (Scene scene, ref Bounds sceneBounds) {
			MeshFilter[] filters = GetComponentsInChildren<MeshFilter> ();
			if (filters == null)
				return false;

			List<Mesh> meshes = new List<Mesh> ();
			List<Instance> instances = new List<Instance> ();
			Dictionary<UnityEngine.Mesh, int> meshToID = new Dictionary<UnityEngine.Mesh, int> ();

			int numMeshes = 0;
			for (int f = 0; f < filters.Length; f++) {

				var uMesh = filters[f].sharedMesh;
				int meshId = -1;
				if (!meshToID.TryGetValue (uMesh, out meshId)) {
					var mesh = Mesh.Create (uMesh, config);
					meshId = numMeshes;
					meshes.Add (mesh);
					meshToID.Add (uMesh, numMeshes++);
				}
				Instance instance = new Instance ();
				instance.meshIndex = meshId;
				instance.storageIdentifier = 0;
				instance.matrix = filters[f].transform.localToWorldMatrix;

				TransformBounds (instance.matrix, meshes[meshId].bounds, ref instance.bounds);
				sceneBounds.Encapsulate (instance.bounds);
				instances.Add (instance);
			}

			scene.meshes = meshes.ToArray ();
			scene.numMeshes = meshes.Count;
			scene.instances = instances.ToArray ();
			scene.numInstances = instances.Count;
			scene.bounds = sceneBounds;

			// #if DEBUG
			// 			DebugExtension.DebugBounds (sceneBounds, 5, false);
			// #endif
			return true;
		}

		[ContextMenu ("Create Default Config")]
		public void CreateDefaultConfig () {
			config = BakeLabConfig.Create ();
		}

		/// <summary>
		/// Bake Signals to vertex attributes
		/// </summary>
		/// <returns>1 - if no errors found</returns>
		[ContextMenu ("Bake")]
		public int Bake () {

			if (config == null) {
				config = BakeLabConfig.Create ();
			}

#if DEBUG
#if UNITY_EDITOR
			UnityEditor.EditorUtility.DisplayProgressBar ("BelfegnarInc.Baking", "Load scene ...", 0.5f);
#endif
			Debug.Log ("Load scene ...");
			Timer totalTimer = new Timer ();
			Timer timer = new Timer ();
			totalTimer.Start ();
			timer.Start ();
#endif

			#region Load scene
			Scene scene = new BelfegnarInc.BakeLab.Scene ();
			Bounds sceneBounds = new Bounds () { min = new Vector3 (float.MaxValue, float.MaxValue, float.MaxValue), max = new Vector3 (float.MinValue, float.MinValue, float.MinValue) };

			if (!LoadScene (scene, ref sceneBounds)) {
#if DEBUG
				Debug.Log ("Failed to load scene, exiting");
#if UNITY_EDITOR
				UnityEditor.EditorUtility.ClearProgressBar ();
#endif
#endif
				return -1;
			}
			#endregion
#if DEBUG
			Debug.Log (timer.Elapsed.TotalMilliseconds);

			#region Print scene stats
			{
				Debug.Log ("Loaded scene: ");
				Debug.Log (scene.numMeshes + " meshes, " + scene.numInstances + " instances");
				int numVertices = 0;
				int numTriangles = 0;
				for (int i = 0; i < scene.numMeshes; ++i) {
					numVertices += scene.meshes[i].numVertices;
					numTriangles += scene.meshes[i].numTriangles;
				}
				Debug.Log ("Uninstanced vertices: " + numVertices);
				Debug.Log ("Uninstanced triangles: " + numTriangles);
			}
			#endregion
#endif

			if (config.flipOrientation) {
				for (int m = 0; m < scene.numMeshes; ++m) {
					Mesh mesh = scene.meshes[m];
					for (int i = 0; i < mesh.numTriangles; i++) {
						Swap (ref mesh.triangles[i * 3 + 0], ref mesh.triangles[i * 3 + 2]);
					}
					if (mesh.normals != null) {
						for (int i = 0; i < mesh.numVertices; i++) {
							mesh.normals[i] *= -1.0f;
						}
					}
				}
			}

#if DEBUG
			Debug.Log ("Minimum samples per face: " + config.minSamplesPerFace);
#if UNITY_EDITOR
			UnityEditor.EditorUtility.DisplayProgressBar ("BelfegnarInc.Baking", "Generate sample points ...", 0.5f);
#endif
			Debug.Log ("Generate sample points ...");
			timer.Reset ();
			timer.Start ();
#endif

			#region Generate samples
			int[] numSamplesPerInstance = new int[scene.numInstances];

#if DEBUG && UNITY_EDITOR
			UnityEditor.EditorUtility.DisplayProgressBar ("BelfegnarInc.Baking", "Distribute samples ...", 0.5f);
#endif
			int totalSamples = Sample.DistributeSamples (scene, config.minSamplesPerFace, config.numSamples, numSamplesPerInstance);

			Samples samples = new Samples ();
			AllocateSamples (samples, totalSamples);

#if DEBUG && UNITY_EDITOR
			UnityEditor.EditorUtility.DisplayProgressBar ("BelfegnarInc.Baking", "Sample instances ...", 0.5f);
#endif
			Sample.SampleInstances (scene, numSamplesPerInstance, config.minSamplesPerFace, samples);
			#endregion

#if DEBUG
			Debug.Log (timer.Elapsed.TotalMilliseconds);

			Debug.Log ("Total samples: " + totalSamples); {
				int sqrtNumRays = (int) (Mathf.Sqrt ((float) (config.numRays)) + .5f);
				Debug.Log ("Rays per sample: " + sqrtNumRays * sqrtNumRays);
				Debug.Log ("Total rays: " + totalSamples * sqrtNumRays * sqrtNumRays);
			}
#endif

			#region Evaluate Signal samples 

#if DEBUG
			Debug.Log ("Compute Signal ...             ");
#if UNITY_EDITOR
			UnityEditor.EditorUtility.DisplayProgressBar ("BelfegnarInc.Baking", "Compute Signal ...", 0.5f);
#endif

			timer.Reset ();
			timer.Start ();
#endif

			float sceneMaxdistance;
			float sceneOffset; {
				float sceneScale = Mathf.Max (Mathf.Max (sceneBounds.max.x - sceneBounds.min.x,
						sceneBounds.max.y - sceneBounds.min.y),
					sceneBounds.max.z - sceneBounds.min.z);
				sceneMaxdistance = sceneScale * config.sceneMaxdistanceScale;
				sceneOffset = sceneScale * config.sceneOffsetScale;
				if (config.sceneOffset != 0) {
					sceneOffset = config.sceneOffset;
				}
				if (config.sceneMaxdistance != 0) {
					sceneMaxdistance = config.sceneMaxdistance;
				}
			}

			float[] floatValues = null;
			List<float[]> floatValuesList = new List<float[]> ();
			int numSignals = 1;

			switch (config.signalType) {
				case SignalType.AO:
					floatValues = new float[totalSamples];
					floatValuesList.Add (floatValues);
					BelfegnarInc.BakeLab.Baking.ComputeAO (scene, samples, config.numRays, sceneOffset, sceneMaxdistance, config.useCpu, floatValues);
					break;

				case SignalType.Radiance:
					floatValues = new float[totalSamples];
					floatValuesList.Add (floatValues);
					BelfegnarInc.BakeLab.Baking.ComputeRadiance (scene, samples, config.numRays, sceneOffset, sceneMaxdistance, config.useCpu, floatValues);
					break;

				case SignalType.RadianceAnimated:
					numSignals = 4;
					Quaternion oldRotation = lightTransform.rotation;
					for (int si = 0; si < numSignals; si++) {
						lightTransform.rotation = Quaternion.Lerp (Quaternion.Euler (fromRotation), Quaternion.Euler (toRotation), si * 1f / (numSignals - 1));
						floatValues = new float[totalSamples];
						floatValuesList.Add (floatValues);
						BelfegnarInc.BakeLab.Baking.ComputeRadiance (scene, samples, config.numRays, sceneOffset, sceneMaxdistance, config.useCpu, floatValues);
					}
					lightTransform.rotation = oldRotation;
					break;

				case SignalType.LightDirTangent:
					numSignals = 3;
					for (int si = 0; si < numSignals; si++) {
						floatValues = new float[totalSamples];
						floatValuesList.Add (floatValues);
						BelfegnarInc.BakeLab.Baking.ComputeTangentSpaceLightDirection (scene, samples, config.numRays, sceneOffset, sceneMaxdistance, config.useCpu, floatValues, si);
					}
					break;

				case SignalType.Irradiance:
					numSignals = 3;
					for (int si = 0; si < numSignals; si++) {
						floatValues = new float[totalSamples];
						floatValuesList.Add (floatValues);
					}
					BelfegnarInc.BakeLab.Baking.ComputeIrradiance (scene, samples, config.packStyle, config.numRays, sceneOffset, sceneMaxdistance, config.useCpu, floatValuesList);
					break;
			}
#if DEBUG
			Debug.Log (timer.Elapsed.TotalMilliseconds);
			Debug.Log ("Map Signal to vertices  ...    ");
			timer.Reset ();
			timer.Start ();
#endif
			List<float[][]> signalsList = new List<float[][]> (numSignals);
			float[][] vertexSignals = null;
			for (int si = 0; si < numSignals; si++) {
				vertexSignals = new float[scene.numInstances][];
				for (int i = 0; i < scene.numInstances; ++i) {
					vertexSignals[i] = new float[scene.meshes[scene.instances[i].meshIndex].numVertices];
				}
				signalsList.Add (vertexSignals);
			}

#if DEBUG && UNITY_EDITOR
			UnityEditor.EditorUtility.DisplayProgressBar ("BelfegnarInc.Baking", "Map Signal to vertices ...", 0);
#endif

			for (int si = 0; si < numSignals; si++) {
				vertexSignals = signalsList[si];
				floatValues = floatValuesList[si];
				Filter.MapSignalToVertices (scene, numSamplesPerInstance, samples, floatValues, config.filterMode, config.regularizationWeight, vertexSignals, parallel);
			}
			#endregion

			//This part aims to fix stitching problems. May be slow though
			if (config.fixStitches) {
				for (int si = 0; si < numSignals; si++) {
					vertexSignals = signalsList[si];
					for (int i = 0; i < vertexSignals.Length; i++) {
						var processedVerts = scene.meshes[scene.instances[i].meshIndex].processedVertices;
						for (int j = 0; j < vertexSignals[i].Length; j++) {
							if (processedVerts.ContainsKey (j)) {
								vertexSignals[i][j] = vertexSignals[i][processedVerts[j]];
							}
						}
					}
				}
			}

			if (!string.IsNullOrEmpty (config.outputFilename)) {
#if DEBUG
				Debug.Log ("Save vertex signals ...");
				timer.Reset ();
				timer.Start ();
#endif
				bool saved = SaveResults (config.outputFilename, scene, signalsList);
#if DEBUG
				Debug.Log (timer.Elapsed.TotalMilliseconds);
				if (saved) {
					Debug.Log ("Saved vertex signals to: " + config.outputFilename);
				} else {
					Debug.Log ("Failed to save vertex signals to: " + config.outputFilename);
				}
#endif
			}

			MeshFilter[] filters = GetComponentsInChildren<MeshFilter> ();
			int fn = 0;
			Color emptyColor = new Color (0, 0, 0, 0);
			float packA = 1, packB = 0;
			if (config.signalType == SignalType.LightDirTangent) {
				packA = 0.5f;
				packB = 0.5f;
			}
			foreach (var f in filters) {
				VertexStreamComponent vsc = f.GetComponent<VertexStreamComponent> ();
				if (vsc == null)
					vsc = f.gameObject.AddComponent<VertexStreamComponent> ();
				Color[] colors = new Color[vertexSignals[fn].Length];
				for (int ic = 0; ic < colors.Length; ic++) {
					colors[ic] = emptyColor;
				}

				for (int si = 0; si < numSignals; si++) {
					var vS = signalsList[si][fn];
					for (int ic = 0; ic < colors.Length; ic++) {
						Color c = colors[ic];
						c[si] = vS[ic] * packA + packB;
						colors[ic] = c;
					}
				}
				vsc.colors = colors;
				fn++;
			}

			signalsList.Clear ();
			//			if(vertexSignals != null)
			//			for (int i = 0; i < scene.numInstances; ++i) {
			//				vertexSignals[i] = null;
			//			}
			vertexSignals = null;
			DestroySamples (samples);
#if DEBUG
			Debug.Log ("Total time: " + totalTimer.Elapsed.TotalMilliseconds + " ms");
#endif

#if DEBUG && UNITY_EDITOR
			UnityEditor.EditorUtility.ClearProgressBar ();
#endif
			System.GC.Collect ();
			return 1;
		}

		/// <summary>
		/// Transforms input bounds with Matrix4x4
		/// </summary>
		/// <param name="matrix">Transformation matrix</param>
		/// <param name="inBounds">Input bounds</param>
		/// <param name="outBounds">Output transformed bounds</param>
		void TransformBounds (Matrix4x4 matrix, Bounds inBounds, ref Bounds outBounds) {
			Vector3 outMin = matrix.MultiplyPoint (inBounds.min);
			Vector3 outMax = matrix.MultiplyPoint (inBounds.max);
			outBounds.SetMinMax (outMin, outMax);
		}

		/// <summary>
		/// Allocates inner arrays in Samples object
		/// </summary>
		/// <param name="samples">Samples object</param>
		/// <param name="n">number of samples</param>
		void AllocateSamples (Samples samples, int n) {
			samples.numSamples = n;
			samples.samplePositions = new Vector3[n];
			samples.sampleNormals = new Vector3[n];
			samples.sampleFaceNormals = new Vector3[n];
			samples.sampleTangents = new Vector4[n];
			samples.sampleInfos = new SampleInfo[n];
			samples.sampleOffset = 0;
		}

		/// <summary>
		/// Clears Samples object
		/// </summary>
		/// <param name="samples">Samples object</param>
		void DestroySamples (Samples samples) {
			samples.samplePositions = null;
			samples.sampleNormals = null;
			samples.sampleFaceNormals = null;
			samples.sampleTangents = null;
			samples.sampleInfos = null;
			samples.numSamples = 0;
			samples.sampleOffset = 0;
		}

		/// <summary>
		/// Saves AO per-vertex information to output file
		/// </summary>
		/// <param name="outputfile">Output file pathname</param>
		/// <param name="s">Scene object</param>
		/// <param name="aoVertex">Per-Instance per-vertex AO values</param>
		/// <returns>True, if no errors found</returns>
		bool SaveResults (string outputfile, Scene scene, List<float[][]> signalsVertex) {
			System.IO.StreamWriter file = new System.IO.StreamWriter (outputfile);

			if (file == null) return false;

			ulong numInstances = (ulong) scene.numInstances;
			ulong numVertices = 0;
			ulong numSignals = (ulong) signalsVertex.Count;

			for (int i = 0; i < scene.numInstances; i++) {
				numVertices += (ulong) scene.meshes[scene.instances[i].meshIndex].numVertices;
			}

			// write header
			file.WriteLine (numInstances);
			file.WriteLine (numVertices);
			file.WriteLine (numSignals);

			// write instances
			ulong vertexOffset = 0;
			for (int i = 0; i < scene.numInstances; i++) {
				ulong identifier = scene.instances[i].storageIdentifier;
				file.WriteLine (identifier);
				file.WriteLine (vertexOffset);
				numVertices = (ulong) scene.meshes[scene.instances[i].meshIndex].numVertices;
				file.WriteLine (numVertices);

				vertexOffset += numVertices;
			}

			// write vertices
			for (int i = 0; i < scene.numInstances; i++) {
				numVertices = (ulong) scene.meshes[scene.instances[i].meshIndex].numVertices;
				for (int si = 0; si < signalsVertex.Count; si++) {
					var sigs = signalsVertex[si][i];
					for (uint v = 0; v < numVertices; v++)
						file.WriteLine (sigs[v]);
				}
			}

			file.Close ();

			return true;
		}
	}

	public class Mesh {
		public bool fixStitches;
		public int numVertices;
		public Vector3[] vertices;
		public Vector3[] normals;
		public Vector4[] tangents;
		public int numTriangles;
		public int[] triangles;
		public Bounds bounds;
		public Dictionary<int, int> processedVertices = new Dictionary<int, int> ();
		public MathNet.Numerics.LinearAlgebra.Double.SparseMatrix regularizationMatrix;

		public static Mesh Create (UnityEngine.Mesh uMesh, BakeLabConfig config) {
			Mesh mesh = null;

			if (!Memory.Meshes.TryGetValue (uMesh, out mesh)) {
				mesh = new Mesh ();
				mesh.numVertices = uMesh.vertexCount;

				List<int> totalIndices = new List<int> ();
				mesh.numTriangles = 0;
				for (int t = 0; t < uMesh.subMeshCount; t++) {
					mesh.numTriangles += (int) uMesh.GetIndexCount (t) / 3;
					totalIndices.AddRange (uMesh.GetIndices (t));
				}

				mesh.vertices = uMesh.vertices;
				mesh.normals = uMesh.normals;
				mesh.tangents = uMesh.tangents;

				//This part aims to fix stitching problems. May be slow though
				if (config.fixStitches) {
					mesh.processedVertices.Clear ();
					var processedVerts = mesh.processedVertices;
					float dist, pEps = config.fixStitchesPositionsDelta, nEps = config.fixStitchesNormalsDelta;
					Vector3 v1, v2;
					for (int i = 0, li = mesh.numVertices; i < li; i++) {
#if DEBUG && UNITY_EDITOR
						if (i % 1000 == 0 && UnityEditor.EditorUtility.DisplayCancelableProgressBar ("BelfegnarInc.BakeLab", "Load scene ...", i * 1f / li))
							return null;
#endif
						v1 = mesh.vertices[i];
						if (!processedVerts.ContainsKey (i))
							for (int j = i + 1; j < li; j++) {
								if (!processedVerts.ContainsKey (j)) {
									v2 = mesh.vertices[j];
									dist = (v1.x - v2.x) * (v1.x - v2.x) + (v1.y - v2.y) * (v1.y - v2.y) + (v1.z - v2.z) * (v1.z - v2.z);
									if (dist <= pEps) {
										v1 = mesh.normals[i];
										v2 = mesh.normals[j];
										dist = (v1.x - v2.x) * (v1.x - v2.x) + (v1.y - v2.y) * (v1.y - v2.y) + (v1.z - v2.z) * (v1.z - v2.z);
										if (dist <= nEps) {
											processedVerts.Add (j, i);
										}
									}
								}
							}
					}
					int index;
					for (int i = 0, li = totalIndices.Count; i < li; i++) {
						index = totalIndices[i];
						if (processedVerts.ContainsKey (index)) {
							totalIndices[i] = processedVerts[index];
						}
					}
				}
				mesh.fixStitches = config.fixStitches;
				mesh.triangles = totalIndices.ToArray ();
				mesh.bounds = uMesh.bounds;
				Memory.Meshes.Add (uMesh, mesh);
			}
			if (mesh != null && mesh.fixStitches != config.fixStitches) {

			}
			return mesh;
		}
	}

	public struct Instance {
		public Matrix4x4 matrix; // 4x4 row major
		public ulong storageIdentifier; // for saving the baked results
		public int meshIndex;
		public Bounds bounds;
	}

	public class Scene {
		public Mesh[] meshes;
		public int numMeshes;
		public Instance[] instances;
		public int numInstances;
		public Bounds bounds;
	}

	public enum VertexFilterMode {
		AREA_BASED = 0,
		LEAST_SQUARES,
		INVALID
	}

	public enum SignalType {
		AO,
		Radiance,
		RadianceAnimated,
		LightDirTangent,
		Irradiance
	}
}