using System.Collections;
using Belfegnar.BakeLab;
using UnityEditor;
using UnityEngine;

[CustomEditor (typeof (Belfegnar.BakeLab.BakeLabMain))]
public class BakeLabMainEditor : Editor {
	public override void OnInspectorGUI () {
		BakeLabMain myTarget = (BakeLabMain) target;

		myTarget.config = EditorGUILayout.ObjectField ("Config", myTarget.config, typeof (BakeLabConfig), true) as BakeLabConfig;
		if (myTarget.config != null) {
			myTarget.parallel = EditorGUILayout.Toggle ("Parallel", myTarget.parallel);
			if (myTarget.config.signalType == SignalType.RadianceAnimated) {
				myTarget.lightTransform = EditorGUILayout.ObjectField ("Light Transform", myTarget.lightTransform, typeof (Transform), true) as Transform;
				myTarget.fromRotation = EditorGUILayout.Vector3Field ("Start Light Rotation", myTarget.fromRotation);
				myTarget.toRotation = EditorGUILayout.Vector3Field ("End Light Rotation", myTarget.toRotation);
			}
			if (GUILayout.Button ("Bake")) {
				myTarget.Bake ();
			}
		}
	}
}