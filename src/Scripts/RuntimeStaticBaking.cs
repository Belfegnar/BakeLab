using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class RuntimeStaticBaking : MonoBehaviour {
    public bool BakeAtAwake = false;
    public float WaitForSeconds = 0f;

    void Awake () {
        if (BakeAtAwake)
            StaticBatchingUtility.Combine (gameObject);
    }

    IEnumerator Start () {
        if (!BakeAtAwake) {
            yield return new WaitForSeconds (WaitForSeconds);
            StaticBatchingUtility.Combine (gameObject);
        }
    }
}