#define BUF_B iChannel1
#define BUF_C iChannel0
#define ORGA  iChannel2
#define BUF_A iChannel3
#define BUMPFACTOR 0.01

#define ORGA_UV_SCALE 0.25+tri(iTime*0.1)*0.25

vec3 get_uv(vec3 p)
{
    float rep =ORGA_UV_SCALE, sca = 0.8;
    vec2 uv = polar(p)*rep;
    //uv.x = mod(uv.x,1.)*0.65+0.25, uv.y = mod(uv.y,1.)*0.65+0.25;
    float vignette = smoothstep(0.9,0.6,length(uv));
    return vec3(uv,vignette);
}

float get_bump(vec3 p)
{
    vec3 uv = get_uv(p);
    vec3 t = texture(ORGA,uv.xy*ORGA_UV_SCALE).xyz*0.8*uv.z;
    return t.x*0.5;
}

vec3 get_bump_norm(vec3 n, vec2 uv)
{
    vec3 uaxis = normalize(cross(vec3(0.0,1.0,0.0), n));
	vec3 vaxis = normalize(cross(uaxis, n));
	mat3 mattanspace = mat3
	(
		uaxis,
		vaxis,
		n
	);
    float delta = -1.0/512.0;
    vec3 a = texture(ORGA, uv + vec2(0.0, 0.0)).xyz;
	float A = a.x;//dot(a,a);
    vec3 b = texture(ORGA, uv + vec2(delta, 0.0)).xyz;
	float B = b.x;//dot(b,b);
    vec3 c = texture(ORGA, uv + vec2(0.0, delta)).xyz;
    float C = c.x;//dot(c,c);


	vec3 norm = normalize(vec3(B - A, C - A, .025));

	return normalize(mattanspace * norm);
}


vec2 drum_hit_planet(vec3 p_front, vec3 p_back)
{
    int block = DRUMS_BLOCK_OFFSET;
    float f, b;
    for(int id = 0; id < 8; id++)
    {
        for(int i = 1; i < 8; i++)
        {
            vec4 data = texelFetch(BUF_A,ivec2(i,block+id),0);
            float env = pow(data.x,0.5);//*Drums_fix_amp[y] ;
            vec3 p  = data.yzw;
            float distf = distance(p,p_front), distb = distance(p,p_back);
            distf = smoothstep(3.,0.5,distf)*env, distb = smoothstep(3.,0.5,distb)*env;
            f += distf, b += distb;
        }
    }
    return vec2(f,b);
}

vec2 bass_hit_planet(vec3 p_front, vec3 p_back)
{
    int block = BASS_BLOCK_OFFSET;
    float f, b;
    for(int id = 0; id < 8; id++)
    {
        for(int i = 1; i < 8; i++)
        {
            vec4 data = texelFetch(BUF_A,ivec2(i,block+id),0);
            float env = pow(data.x,0.5);//*Drums_fix_amp[y] ;
            vec3 p  = data.yzw;
            float distf = distance(p,p_front), distb = distance(p,p_back);
            distf = smoothstep(3.,0.5,distf)*env, distb = smoothstep(3.,0.5,distb)*env;
            f += distf, b += distb;
        }
    }
    return vec2(f,b);
}



void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = (2.0*(fragCoord)-iResolution.xy)/iResolution.y;
    
    vec4 bufb = texture(BUF_B,fragCoord/iResolution.xy);
    vec3 col = bufb.xyz, glow = FLOWER_COL_CENTER*bufb.w;
    float fgl = bufb.w ;
    vec3 ro = false ? vec3(cos(iTime),0.5,sin(iTime))*10. : texelFetch(BUF_A,RO_COO,0).xyz;
    vec3 lookat = LOOKAT;
    mat3 cam = camera(ro, lookat, 0.);
    vec3 rd = cam * normalize(vec3(uv,1));
    float rad = STAR_RAD;
    vec2 sph = asphere(ro,rd,STAR_RAD);
    float b = get_bump(ro+rd*sph.x);
    int song_sect = int(texelFetch(BUF_C, ivec2(0),0).x);
    float front=  sph.x, back = sph.y;
    vec3 p_front = ro+rd*front, p_back = ro+rd*back;
    vec3 uv_fr = get_uv(ro+rd*front), uv_bk = get_uv(ro+rd*back);
    vec3 txf = texture(ORGA,uv_fr.xy).xyz*1.0*uv_fr.z;
    vec3 txb = texture(ORGA,uv_bk.xy).xyz*1.0*uv_bk.z;
    vec4 dru = texture(BUF_C,fragCoord/iResolution.xy);
    float d_ctr = sphDensity(ro,rd,vec3(0),10.,100.);
    
    fgl *= clamp(bufb.w*3.,0.,1.);
    if(sph.x < MAX_DIST && sph.x > 0.2)
    {
        //glow*= 0.8;
        fgl = smoothstep(0.2,2.9,fgl*20.)*0.5;//pow(clamp(fgl*0.5,0.,1.),2.)*0.5;
        vec3 norm_front = normalize(ro+rd*sph.x), norm_back = normalize(ro+rd*sph.y);
        // bump mapping
	    vec3 surf_norm_front = get_bump_norm(norm_front,uv_fr.xy), 
             surf_norm_back = get_bump_norm(norm_back , uv_bk.xy);
        
        if((song_sect > 3 && song_sect < 6) || song_sect > 6)
        {
            float fresnel = clamp(1. - dot(surf_norm_front ,-rd), 0.,1.);
            fresnel += clamp(1. - dot(surf_norm_back ,-rd), 0.,1.);
        
            col += SunCol*pow(fresnel,0.1)*0.01;
        }
     
        col = max(blinn_phong(ro+rd*front, rd, vec3(0,1,0),surf_norm_front, SunCol)*fgl,col);
        col = max(blinn_phong(ro+rd*back , rd, vec3(0,1,0),surf_norm_back, SunCol )*fgl,col);
        vec2 drum_hit = drum_hit_planet(p_front,p_back);
        vec2 bass_hit = bass_hit_planet(p_front,p_back);
        col += (txf*drum_hit.x+txb*drum_hit.y)*SunCol*5.;
        col += (txf*bass_hit.x+txb*bass_hit.y)*vec3(0.9647, 0.0745, 0.0745)*5.;
        //dru.xyz *= smoothstep(0.9,0.5,d_ctr);
        //col = vec3(fgl);
    }


    #if 0
    vec2 o = vec2(0.0005,0);
    col += texture(BUF_B,fragCoord/iResolution.xy+o).xyz;
    col += texture(BUF_B,fragCoord/iResolution.xy+o.yx).xyz;
    col += texture(BUF_B,fragCoord/iResolution.xy-o).xyz;
    col += texture(BUF_B,fragCoord/iResolution.xy-o.yx).xyz;
    col /= 4.;
    #endif
    col = pow(col+clamp(glow,0.,1.),vec3(.4545));
    
    #if 0
    //fragColor = vec4(getSongSection(MIDI) == 2 ? 1 : 0);
    
    fragColor = true ?  texelFetch(BUF_A,ivec2(fragCoord*vec2(0.02,0.51)+vec2(0,0)),0) 
                    : texelFetch(BUF_C,ivec2(fragCoord*vec2(1)),0).xxxx;
    
    #else
    
    //col = dot(col,col) > 0.01 ? min(vec3(glow),col): col;
    //glow *= dot(col,col) > 0.1 ? 0. : 2.;
    fragColor = vec4(col+dru.xyz,0.);
    //fragColor = vec4(col,0.);
    #endif
}