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

vec4 drawDrums(vec3 ro, vec3 rd)
{
    float glow;
    vec3 col; 
    vec3 hit_point = vec3(0);
    int block = DRUMS_BLOCK_OFFSET+DRUMS_BLOCK.z;
    for(int id = 0; id < 8; id++)
    {
        float env = pow(tri(texelFetch(BUF_A,ivec2(id,DRUMS_ENV_ROW),0).x),3.) ;
        vec3 p  = texelFetch(BUF_A,ivec2(0,id+block),0).xyz;
        
        //pp.xy *=1.;
        glow += integrateLightFullView(ro-p,rd,5.1*env,0.1);// (100.8*pow(env,0.5))/(0.1+pa*pa*50.);
        vec3 ha = hash31(float(id));
        col += glow*ha;
    }
    return vec4(col,glow);
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
    int song_sect = getSongSection(MIDI);
    bool drums_on = is_drums_on(song_sect);
    vec4 dru = drums_on ? drawDrums(ro,rd) : vec4(0);
    vec4 feedb = texture(FEEDBACK,fragCoord/iResolution.xy);
    vec4 col = dru +feedb*0.5;
    fragColor = col;
    

}