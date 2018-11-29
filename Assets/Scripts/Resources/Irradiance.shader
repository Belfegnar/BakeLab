// Upgrade NOTE: replaced '_World2Object' with 'unity_WorldToObject'

// Upgrade NOTE: replaced 'mul(UNITY_MATRIX_MVP,*)' with 'UnityObjectToClipPos(*)'

Shader "BelfegnarInc/Baking/Irradiance" {
	Properties {
		_Color("Color", Color) = (1,1,1,1)
		_MainTex("Albedo", 2D) = "white" {}
		_BumpMap ("Normal (XYZ*0.5+0.5)", 2D) = "bump" {}
	}
	SubShader {
		Tags { "RenderType"="Opaque" }
		ZWrite On 
		Cull Back
		LOD 100
		Pass{
			Tags {"LightMode"="ForwardBase"}
			CGPROGRAM
			#pragma vertex vert
			#pragma fragment frag
			#include "UnityCG.cginc"
			#include "UnityLightingCommon.cginc" // for _LightColor0
	
			struct appdata
			{
				float4 vertex : POSITION;
				float3 normal : NORMAL;
				float2 texcoord : TEXCOORD0;
				float4 occl : COLOR;
			};
			
			struct v2f
			{
				float4 vertex : SV_POSITION;
				float2 texcoord : TEXCOORD0;
				float3 color0 : TEXCOORD1;
			};
			
			fixed4 _Color;
			sampler2D _MainTex;
			sampler2D _BumpMap;
			float4 _MainTex_ST;
			
			v2f vert (appdata v)
			{
				v2f o;
				o.vertex = UnityObjectToClipPos(v.vertex);
				o.texcoord = TRANSFORM_TEX(v.texcoord, _MainTex);
				o.color0 = pow(max(0, v.occl.rgb),1.0/2.2);
				// Uncomment for direct lighting
				// half3 worldNormal = UnityObjectToWorldNormal(v.normal);
				// half nl = max(0.05, dot(worldNormal, _WorldSpaceLightPos0.xyz));
                // o.color0 += nl * _LightColor0.rgb*_LightColor0.a;
				return o;
			}
	
			fixed4 frag (v2f i) : Color
			{
				fixed4 albedo = fixed4(1,1,1,1);// tex2D(_MainTex, i.texcoord)*_Color;
				return fixed4(i.color0*albedo.rgb, 1);
			}
			ENDCG
		}
	} 
	//FallBack "Diffuse"
}