using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace BelfegnarInc.BakeLab {
	/// <summary>
	/// Contains configuration settings.
	/// </summary>
	[CreateAssetMenu (fileName = "BakeLabConfig", menuName = "BakeLab/Config", order = 1)]
	public class BakeLabConfig : ScriptableObject {

		const int NUM_RAYS = 64;
		const int SAMPLES_PER_FACE = 3;
		const float GROUND_SCALE = 100.0f;
		const float GROUND_OFFSET = 0.03f;
		const float SCENE_OFFSET_SCALE = 0.001f;
		const float SCENE_MAXDISTANCE_SCALE = 1.1f;
		const float REGULARIZATION_WEIGHT = 0.5f;

		public SignalType signalType = SignalType.AO;
		[HideInInspector]
		public BakeLab.Baking.PackStyle packStyle = BakeLab.Baking.PackStyle.Dir3;
		public int numSamples = 0; // default means determine from mesh
		public int minSamplesPerFace = SAMPLES_PER_FACE;
		public int numRays = NUM_RAYS;
		public VertexFilterMode filterMode = VertexFilterMode.LEAST_SQUARES;
		public float regularizationWeight = REGULARIZATION_WEIGHT;
		public float sceneOffsetScale = SCENE_OFFSET_SCALE;
		public float sceneMaxdistanceScale = SCENE_MAXDISTANCE_SCALE;
		public float sceneMaxdistance = 0;
		public float sceneOffset = 0; // must default to 0
		public bool useCpu = false;
		public bool flipOrientation = false;
		public string outputFilename;
		public bool fixStitches = true;
		public float fixStitchesPositionsDelta = float.Epsilon;
		public float fixStitchesNormalsDelta = float.Epsilon;

		[ContextMenu ("Set Defaults")]
		public void SetDefaults () {

			numSamples = 0; // default means determine from mesh
			minSamplesPerFace = SAMPLES_PER_FACE;
			numRays = NUM_RAYS;
			sceneOffsetScale = SCENE_OFFSET_SCALE;
			sceneMaxdistanceScale = SCENE_MAXDISTANCE_SCALE;
			sceneOffset = 0; // must default to 0
			sceneMaxdistance = 0;
			useCpu = false;
			flipOrientation = false;
			filterMode = VertexFilterMode.LEAST_SQUARES;
			regularizationWeight = REGULARIZATION_WEIGHT;
		}

		public static BakeLabConfig Create () {
			BakeLabConfig result = ScriptableObject.CreateInstance (typeof (BakeLabConfig)) as BakeLabConfig;
			result.SetDefaults ();
			return result;
		}
	}
}