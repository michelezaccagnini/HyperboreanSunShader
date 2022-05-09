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
    //vec4 t;
    vec3 hit_point = vec3(0);
    int block = DRUMS_BLOCK_OFFSET;

    for(int id = 0; id < 8; id++)
    {   
        // 7 - kick, 4 - clave a, 2 - clave b, 6 - snare, 1 -hats,
        //if(id != 7 ) continue;
        //I am not using some of the drums hits
        //int id = Relevant_drums[y];
        //if(id !=6) continue;
        for(int i = 1; i < 8; i++)
        {
            vec4 data = texelFetch(BUF_A,ivec2(i,block+id),0);
            float env = pow(data.x*1.,0.65);//*Drums_fix_amp[y];
            vec3 p  = data.yzw;
            glow += integrateLightFullView(ro-p,rd,5.1*env,0.051)*env;// (100.8*pow(env,0.5))/(0.1+pa*pa*50.);
            vec3 ha = hash31(float(id))*SunCol*1.;
            col += glow*ha;
        }
       
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
    vec4 col = dru*0.8 +feedb*0.5;
    fragColor = col;
    

}