//Bass buffer red lights bouncing on planet
#define FEEDBACK iChannel0
#define BUF_A iChannel1
#define MIDI iChannel2

vec4 drawBass(vec3 ro, vec3 rd)
{
    float glow;
    vec3 col; 
    //vec4 t;
    vec3 hit_point = vec3(0);
    int block = BASS_BLOCK_OFFSET;

    for(int id = 0; id < 8; id++)
    {   
        
        //I am not using some of the drums hits
        //if(id !=6) continue;
        for(int i = 1; i < 8; i++)
        {
            vec4 data = texelFetch(BUF_A,ivec2(i,block+id),0);
            float env = pow(data.x*1.,0.65);//*Drums_fix_amp[y];
            vec3 p  = data.yzw;
            glow += integrateLightFullView(ro-p,rd,10.1*env,0.041)*env;// (100.8*pow(env,0.5))/(0.1+pa*pa*50.);
            vec3 ha = hash31(float(id*5.))*vec3(0.9647, 0.051, 0.051);
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
    bool bass_on = is_bass_on(song_sect);
    vec4 bas = bass_on ? drawBass(ro,rd) : vec4(0);
    vec4 feedb = texture(FEEDBACK,fragCoord/iResolution.xy);
    vec4 col = bas*0.6 +feedb*0.3;
    fragColor = col;
    

}