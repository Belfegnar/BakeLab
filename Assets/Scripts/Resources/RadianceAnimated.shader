// Upgrade NOTE: replaced '_World2Object' with 'unity_WorldToObject'

// Upgrade NOTE: replaced 'mul(UNITY_MATRIX_MVP,*)' with 'UnityObjectToClipPos(*)'

Shader "BelfegnarInc/Baking/RadianceAnimated" {
	Properties {
		_MainTex ("Base (RGB)", 2D) = "white" {}
		_BumpMap ("Normal (XYZ*0.5+0.5)", 2D) = "bump" {}
	}
	SubShader {
		Tags { "RenderType"="Opaque" }
		ZWrite On 
		Cull Back
		LOD 100
		Pass{
			CGPROGRAM
			#pragma vertex vert
			#pragma fragment frag
			#include "UnityCG.cginc"
	
			struct appdata
			{
				float4 vertex : POSITION;
				float2 texcoord : TEXCOORD0;
				float4 occl : COLOR;
			};
			
			struct v2f
			{
				float4 vertex : SV_POSITION;
				float2 texcoord : TEXCOORD0;
				float3 color : COLOR;
			};
			
			sampler2D _MainTex;
			sampler2D _BumpMap;
			float4 _MainTex_ST;
			uniform float4 _Animation;
			
			v2f vert (appdata v)
			{
				v2f o;
				o.vertex = UnityObjectToClipPos(v.vertex);
				o.texcoord = TRANSFORM_TEX(v.texcoord, _MainTex);
				o.color = (dot(v.occl, _Animation)).xxx;
				return o;
			}
	
			fixed4 frag (v2f i) : Color
			{
				return fixed4(i.color, 1);
			}
			ENDCG
		}
	} 
	//FallBack "Diffuse"
}