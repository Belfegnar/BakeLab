using System.Collections;
using System.Collections.Generic;
using UnityEngine;

[ExecuteInEditMode]
public class RadianceAnimator : MonoBehaviour {

	public float Speed = 0.2f;
	[SerializeField] int Frame = 0;
	[SerializeField] int NextFrame = 0;
	[SerializeField] float Diff;
	[SerializeField] Vector4 AnimationValue;

	void OnEnable () {
		Shader.SetGlobalVector ("_Animation", new Vector4 (1, 0, 0, 0));
	}

	void OnDisable () {
		Shader.SetGlobalVector ("_Animation", new Vector4 (1, 0, 0, 0));
	}

	void Update () {
		var count = Time.time * Speed;
		var countInt = (int) count;
		Diff = count - countInt;
		Frame = (int) (Time.time * Speed) % 4;
		NextFrame = Frame >= 3 ? 0 : (Frame + 1);
		int index = -1;
		while (index++ < 3) {
			float v = 0f;
			if (index == Frame)
				v = 1f - Diff;
			else if (index == NextFrame)
				v = Diff;
			AnimationValue[index] = v;
		}
		Shader.SetGlobalVector ("_Animation", AnimationValue);
	}
}