using System.Collections;
using BelfegnarInc.BakeLab;
using UnityEditor;
using UnityEngine;

[CustomEditor (typeof (BelfegnarInc.BakeLab.BakeLabMain))]
public class BakeLabMainEditor : Editor {
	public override void OnInspectorGUI () {
		BakeLabMain myTarget = (BakeLabMain) target;

		myTarget.config = EditorGUILayout.ObjectField ("Config", myTarget.config, typeof (BakeLabConfig)) as BakeLabConfig;
		if (myTarget.config != null) {
			myTarget.parallel = EditorGUILayout.Toggle ("Parallel", myTarget.parallel);
			if (myTarget.config.signalType == SignalType.RadianceAnimated) {
				myTarget.lightTransform = EditorGUILayout.ObjectField ("Light Transform", myTarget.lightTransform, typeof (Transform)) as Transform;
				myTarget.fromRotation = EditorGUILayout.Vector3Field ("Start Light Rotation", myTarget.fromRotation);
				myTarget.toRotation = EditorGUILayout.Vector3Field ("End Light Rotation", myTarget.toRotation);
			}
			if (GUILayout.Button ("Bake")) {
				myTarget.Bake ();
			}
		}
	}
}