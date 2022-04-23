#define BUF_B iChannel1
#define BUF_A iChannel0
#define BUF_C iChannel2
#define MIDI iChannel3
void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = (2.0*(fragCoord)-iResolution.xy)/iResolution.y;
    
    vec3 col = texture(BUF_B,fragCoord/iResolution.xy).xyz;//max(texture(BUF_C,fragCoord/iResolution.xy).xyz,texture(BUF_B,fragCoord/iResolution.xy).xyz); 
    vec3 ro = false ? vec3(cos(iTime),0.5,sin(iTime))*20. : texelFetch(BUF_A,RO_COO,0).xyz;
    vec3 lookat = vec3(0);
    mat3 cam = camera(ro, lookat, 0.);
    vec3 rd = cam * normalize(vec3(uv,1));
    float rad = STAR_RAD;
    vec2 sph = asphere(ro,rd,STAR_RAD);
    float front=  sph.x, back = sph.y;
    vec2 sph_uv1 = polar(normalize(ro+rd*front));
    vec2 sph_uv2 = polar(ro+rd*back);
    if(sph.x < 80. && sph.x > 0.2)
    {
        //col += max(texture(BUF_B,sph_uv1-vec2(0.5,0)).xyz*0.8,vec3(0));
        col += texture(BUF_B,sph_uv2).xyz*0.8;
        //col = vec3(sph_uv1.xxx);
    }


    /*
    vec2 o = vec2(0.0015,0);
    col += texture(BUF_B,fragCoord/iResolution.xy+o).xyz;
    col += texture(BUF_B,fragCoord/iResolution.xy+o.yx).xyz;
    col += texture(BUF_B,fragCoord/iResolution.xy-o).xyz;
    col += texture(BUF_B,fragCoord/iResolution.xy-o.yx).xyz;
    col /= 4.;
    */
    #if 0
    //fragColor = vec4(getSongSection(MIDI) == 2 ? 1 : 0);
    
    fragColor = false ? texelFetch(BUF_C,ivec2(fragCoord*vec2(1)),0) : 
                    texelFetch(BUF_A,ivec2(fragCoord*vec2(0.03,0.251)),0).xxxx  ;
    
    #else
    fragColor = vec4(col,0.);
    #endif
}