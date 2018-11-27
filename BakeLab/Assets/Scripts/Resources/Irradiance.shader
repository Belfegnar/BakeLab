// Upgrade NOTE: replaced '_World2Object' with 'unity_WorldToObject'

// Upgrade NOTE: replaced 'mul(UNITY_MATRIX_MVP,*)' with 'UnityObjectToClipPos(*)'

Shader "BelfegnarInc/Baking/Irradiance" {
	Properties {
		_MainTex ("Base (RGB)", 2D) = "white" {}
		_NormalTex ("Normal (XYZ*0.5+0.5)", 2D) = "gray" {}
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
				//float3 normal : NORMAL;
				float2 texcoord : TEXCOORD0;
				float4 occl : COLOR;
			};
			
			struct v2f
			{
				float4 vertex : SV_POSITION;
				float2 texcoord : TEXCOORD0;
				float3 color0 : TEXCOORD1;
			};
			
			sampler2D _MainTex;
			sampler2D _NormalTex;
			float4 _MainTex_ST;
			
			v2f vert (appdata v)
			{
				v2f o;
				o.vertex = UnityObjectToClipPos(v.vertex);
				o.texcoord = TRANSFORM_TEX(v.texcoord, _MainTex);
				o.color0 = v.occl.rgb;
				return o;
			}
	
			fixed4 frag (v2f i) : Color
			{
				return fixed4(pow(i.color0, 1.0/2.2), 1);
			}
			ENDCG
		}
	} 
	//FallBack "Diffuse"
}