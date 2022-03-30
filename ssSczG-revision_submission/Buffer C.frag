//MIDI data parsing
#define FEEDBACK iChannel0
#define BUF_A iChannel1
#define MIDI iChannel2

vec3 drawTrails(vec3 ro, vec3 rd, int num_trails, int tex_coo, int env_row, sampler2D text)
{
    vec3 lig = vec3(0);
    for(int id = 0; id < num_trails; id++)
    {
        vec3 pos = texelFetch(text,ivec2(0,tex_coo+id),0).xyz;
        float env = pow(texelFetch(text,ivec2(id, env_row),0).x,1.5);
        float l = smoothstep(0.1,0.99,integrateLightFullView(ro-pos,rd,15.01,0.041));
        lig += max(l*vec3(1)*env, vec3(0));
    }
    return  lig;
    
}

vec3 drawDots(vec2 uv, int num_trails, int tex_coo, int env_row, sampler2D text)
{
    vec3 res = vec3(0);
    for(int id = 0; id < num_trails; id++)
    {
        vec2 pos = texelFetch(text,ivec2(0,tex_coo+id),0).xy;
        float env = pow(texelFetch(text,ivec2(id, env_row),0).x,1.);
        float l = smoothstep(0.08,0.07,length(uv-pos));
        res += max(l*vec3(1)*env, vec3(0));
    }
    return  res;
    
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
     // Normalized pixel coordinates (from 0 to 1)
    vec2 uv = (fragCoord-0.5*iResolution.xy)/iResolution.y*2.;
    vec2 uv2 = fragCoord/iResolution.y;
    vec3 ro = texelFetch(BUF_A,RO_COO,0).xyz;
    vec3 lookat = vec3(0);
    mat3 cam = camera(ro, lookat, 0.);
    vec3 rd = cam * normalize(vec3(uv,1));
    vec3 lig = drawDots(uv2,8, PERHI_BLOCK.w, PERHI_ENV_ROW, BUF_A);
    vec3 feedb = texture(FEEDBACK,fragCoord/iResolution.xy).xyz;
    vec3 col = lig +feedb*0.1;
    fragColor.xyz = col;
    

}