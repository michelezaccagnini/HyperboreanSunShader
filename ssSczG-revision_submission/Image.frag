#define BUF_B iChannel1
#define BUF_A iChannel0
#define BUF_C iChannel2
void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec3 col = texture(BUF_B,fragCoord/iResolution.xy).xyz;//max(texture(BUF_C,fragCoord/iResolution.xy).xyz,texture(BUF_B,fragCoord/iResolution.xy).xyz); 
    #if 0
    fragColor = true ? texelFetch(BUF_A,ivec2(fragCoord*vec2(0.01,0.2)),0) : 
                    abs(texelFetch(BUF_A,ivec2(.038*fragCoord),0)).x > 0. ? vec4(1) : vec4(0)  ;
    #else
    fragColor = vec4(col,0.);
    #endif
}