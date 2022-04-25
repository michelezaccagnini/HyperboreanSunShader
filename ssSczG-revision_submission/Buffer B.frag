#define FEEDBACK iChannel0
#define BUF_A iChannel1
#define MIDI iChannel2
#define TEXTURE iChannel3
#define ZERO (min(iFrame,0))
#define AA 1


/*
musical data explanation:
objects move back and forth (bouncing motion) on 3 dimensions
each one-dimensional bouncing motion can have a different period
at a given moment their position is "measured" (snapshot of their positions) 
and their parameters are returned determining the position in time, 
velocity and pitch: while their movement is continuous, 
their sound parameters are output only at certain points in time (not randomly)

each dimension is mapped to a parameter: 
1. rhythmic position (where in bar)
2. velocity
3. pitch
*/
#define MAX_DIST 100.
struct HitInfo
{
    float dist;
    //id is element id  and sub element id (e.g. rope id)
    ivec2 id;
    //returns rope/sphere connection interpolation index
    float smin;
    vec3 pos;
    vec3 surf;
    vec3 nor;
    vec2 uv;
    vec2 uv_transorm;
    vec3 col;
    float env;
};

float glow;

HitInfo map(vec3 p)
{
    HitInfo res;
    res.dist = MAX_DIST;
    int song_sect = getSongSection(MIDI);
    bool flower_on = is_flower_on(song_sect);
    bool pong_on = is_pong_on(song_sect);
    bool perhi_on = is_perhi_on(song_sect);
    bool drums_on = is_drums_on(song_sect);
    bool bass_on = is_bass_on(song_sect);
    vec3 hit_point = vec3(0);
    //flower
    if(flower_on)
    {
        vec3 hit_point = vec3(0),
             hit_point2 = vec3(0);
        for(int id = 0; id < 5; id++)
        {
            float env = texelFetch(BUF_A,ivec2(id,FLOWER_ENV_ROW),0).x;
            vec3 displ = flower_sect_displ(song_sect);
            vec3 sp_pos = displ;
            vec3 rop = rope_flower1(p-displ,id+FLOWER_BLOCK.z  ,BUF_A,env,hit_point);
            //vec3 rope_attach = rope_flower2(p-displ,id+FLOWER_BLOCK.z*3,BUF_A,hit_point2); 
            //float joint = length(p-texelFetch(BUF_A,ivec2(0,id+FLOWER_BLOCK.z*3),0).xyz)-0.1;
            //vec3 rop = opU(rope_center,rope_attach);
            //bool is_first = rope_center.x < rope_attach.x; 
            float c_sphere = length(p-sp_pos)-0.3;
            hit_point = hit_point;//is_first ? hit_point : hit_point2;
            float h = 0.;
            rop.x = smin(c_sphere,rop.x,0.7,h);
            float glow_int = smoothstep(0.051,0.,abs(rop.z*0.5-(1.-pow(env,0.5))))*0.1*env;
            glow += (0.1/(0.1+rop.x*rop.x*rop.x))*glow_int;
            if(rop.x < res.dist)
            {
                //id += is_first ? 0 : 5;
                vec2 tuv = rop.yz*0.1;
                tuv.y += env;
                vec3 txt = texture(TEXTURE,tuv).xyz;
                float bump = clamp(dot(txt,txt),0.,1.)*0.051;
                res.dist = rop.x-bump, res.id = ivec2(1,id), res.uv = rop.yz,
                res.uv_transorm = tuv, 
                res.pos = p, res.surf = hit_point,
                res.env = env, 
                res.smin = h;           
            }
        }
    }
    if(pong_on)
    {
        vec3 hit_point = vec3(0);
        for(int id = 0; id < 8; id++)
        {
            float env = texelFetch(BUF_A,ivec2(id,PONG_ENV_ROW),0).x;
            //float sph = length(p - pos) - 0.4;
            vec3 rop  = ropePong(p, id+PONG_BLOCK_OFFSET+PONG_BLOCK.z,BUF_A,env, hit_point);
            float glow_int = smoothstep(0.051,0.,abs(rop.z*0.5-env))*0.5*env;
            glow += (0.1/(0.1+rop.x*rop.x))*glow_int*FLOWER_COL_CENTER;
            //rop.x = min(rop.x,sph);
            if(rop.x < res.dist)
            {
                vec2 tuv = rop.yz*0.1;
                tuv.y += env;
                vec3 txt = texture(TEXTURE,tuv).xyz;
                float bump = clamp(dot(txt,txt),0.,1.)*0.01;
                res.dist = rop.x-bump, res.id = ivec2(2,id), 
                res.uv_transorm = tuv, res.uv = rop.yz, res.pos = p, res.surf = hit_point, 
                res.env = texelFetch(BUF_A,ivec2(id,PONG_ENV_ROW),0).x,
                res.surf = hit_point, res.smin = 0.;
            }
        }
    }
    if(perhi_on )
    {
        ivec2 o_sect = any(equal(ivec2(song_sect),ivec2(2,4))) ? ivec2(0 ,4) : ivec2(0,8) ; 
        vec3 hit_point = vec3(0);
        for(int id = o_sect.x; id < o_sect.y; id++)
        {
            float env = texelFetch(BUF_A,ivec2(id,PERHI_ENV_ROW),0).x;
            vec3 displ = perhi_sect_displ(song_sect);
            vec3 sp_pos = displ ;
            float center_sph = length(p-sp_pos) -0.7;
            vec3 rop  = ropePerhi(p-displ, id+PERHI_BLOCK_OFFSET+PERHI_BLOCK.z,BUF_A, env, hit_point);
            //rop.x = min(rop.x,sph); 
            float h = 0.;
            float glow_int = smoothstep(0.051,0.,abs(rop.z*0.5-pow(env,0.5)))*0.5*pow(env,.85);
            glow += (0.1/(0.1+rop.x*rop.x))*glow_int*PERHI_COL_CENTER;
            rop.x = smin(rop.x,center_sph,0.6,h);
            if(rop.x < res.dist)
            {
                vec2 tuv = rop.yz*0.1;
                tuv.y += env;
                vec3 txt = texture(TEXTURE,tuv).xyz;
                float bump = clamp(dot(txt,txt),0.,1.)*0.1;
                res.dist = rop.x-bump, res.id = ivec2(3,id), res.uv_transorm = tuv, 
                res.uv = rop.yz, res.pos = p, 
                res.surf = hit_point, 
                res.env = texelFetch(BUF_A,ivec2(id,PERHI_ENV_ROW),0).x,
                res.surf = hit_point,
                res.smin = 1.-h;
            }
        }
    }
    if(drums_on)
    {
        vec3 hit_point = vec3(0);
        int block = DRUMS_BLOCK_OFFSET+DRUMS_BLOCK.z;
        for(int id = 0; id < 8; id++)
        {
            float env = texelFetch(BUF_A,ivec2(id,DRUMS_ENV_ROW),0).x;
            vec3 rop  = ropeDrums(p, id+block,BUF_A, env, hit_point);
            
            //rop.x = min(rop.x,sph);
            if(rop.x < res.dist)
            {
                vec2 tuv = rop.yz*0.1;
                tuv.y += env;
                vec3 txt = texture(TEXTURE,tuv).xyz;
                float bump = clamp(dot(txt,txt),0.,1.)*0.1;
                res.dist = rop.x, res.id = ivec2(4,id), res.uv_transorm = tuv,
                res.uv = rop.yz, 
                res.pos = p, res.surf = hit_point, 
                res.env = texelFetch(BUF_A,ivec2(id,DRUMS_ENV_ROW),0).x;
                res.surf = hit_point, 
                res.smin = 0.;
            }
        }
    }
    if(bass_on)
    {
        vec3 hit_point = vec3(0);
        for(int id = 0; id < 8; id++)
        {
            float env = texelFetch(BUF_A,ivec2(id,BASS_ENV_ROW),0).x;
            vec3 displ = song_sect > 5 ?vec3(0,0, LAST_SECT_DISPLACE) : vec3(0);
            vec3 rop  = ropeDrums(p-displ, id+BASS_BLOCK_OFFSET+BASS_BLOCK.z,BUF_A, env, hit_point);
            //rop.x = min(rop.x,sph);
            if(rop.x < res.dist)
            {
                res.dist = rop.x, res.id = ivec2(4,id), res.uv = rop.yz, 
                res.pos = p, res.surf = hit_point, 
                res.env = texelFetch(BUF_A,ivec2(id,BASS_ENV_ROW),0).x;
                res.surf = hit_point, res.smin = 0.;
            }
        }
    }
    //res.dist *= 0.7;
    return res;
}

HitInfo intersect(vec3 ro, vec3 rd)
{
    HitInfo res;
    float min_dist = 0.01;
    float max_dist = MAX_DIST;
    float d = min_dist;
    res.id = ivec2(-1);
    for(int i = 0; i < 80; i++)
    {
        vec3 p = ro + rd*d;
        res = map(p);
        if(d > max_dist || abs(res.dist) < min_dist) break;
        d += res.dist;
    }
    res.dist = d;
    return res;
}

float intersectSphere(vec3 ro, vec3 rd, vec3 pos, float rad)
{
    float min_dist = 0.01;
    float max_dist = MAX_DIST;
    float d = min_dist;
    for(int i = 0; i < 80; i++)
    {
        vec3 p = ro + rd*d;
        float s = length(p-pos)-rad;
        if(d > max_dist || s < min_dist) break;
        d += s;
    }
    return d;
}


vec3 normal( in vec3 pos)
{
    const float ep = 5.e-4;
    vec2 e = vec2(1.,0.);
    float s = map(pos).dist;
    return normalize( vec3(map( pos + e.xyy*ep).dist, 
                           map( pos + e.yxy*ep).dist, 
                           map( pos + e.yyx*ep).dist) - s);
}

vec3 calcLight(HitInfo hit, vec3 rd, vec3 lig_pos)
{ 
    vec3 diff = max(0.,dot(hit.nor, lig_pos)*0.5+0.5)*vec3(0.2,1.12,1.26);
    //reflect
    vec3 refl = reflect(lig_pos,hit.pos);
    vec3 spec_lig = vec3(clamp(dot(rd, refl),0.,1.));
    int song_sect =getSongSection(MIDI);
    vec3 sun_pos = vec3(0), sun_col = vec3(0);
    if(any(equal(ivec3(hit.id.x),ivec3(1,2,4))))
    {
        sun_pos = flower_sect_displ(song_sect);
        sun_col = FLOWER_COL_CENTER;
    }
    else
    {
        sun_pos = perhi_sect_displ(song_sect);
        sun_col = PERHI_COL_CENTER;
    }
    float light_distance = smoothstep(.1,2., length(sun_pos- hit.pos)*0.2);
    float falloff = .05;
    vec3 ha = hash31(float(hit.id.x)+float(hit.id.y));
    hit.col  = max(ha*1.2+sun_col,vec3(1));
    float light_intensity = clamp(falloff/pow(light_distance,1.8),0,1);
    float fres = clamp(1. - dot(hit.nor, -rd),0.,1.);
    vec3 col = fres*0.2*diff*hit.col*light_intensity;
    float stripe = 0.;
    if(hit.id.x == 3) hit.env = 1.- hit.env;
    stripe = pow(smoothstep(0.1,0.0,abs(pow(hit.env,0.5)-(1.-hit.uv.y))),0.5);
    col *=vec3(stripe *3.+0.8);//diff*pow(hit.env,1.5)*1.+0.1;//max(col, vec3(0));
    vec3 txt = texture(TEXTURE,hit.uv_transorm).xyz;
    float bump = clamp(dot(txt,txt),0.,1.);
    col += txt*light_intensity*0.5;
    col += bump*PERHI_COL_CENTER*0.2;
    // color center spheres for FLOWER and PERHI
    col =  mix(col,sun_col,hit.smin);
    if(hit.id.x == 1 && hit.id.y > 4) col *= 0.2;
    return col;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec3 tot = vec3(0);
#if AA > 1
    for( int m=ZERO; m<AA; m++ )
    for( int n=ZERO; n<AA; n++ )
    {
        // pixel coordinates
        vec2 o = vec2(float(m),float(n)) / float(AA)  ;
        vec2 uv = (2.0*(fragCoord+o)-iResolution.xy)/iResolution.y;
#else
        vec2 uv = (2.0*(fragCoord)-iResolution.xy)/iResolution.y;
#endif
        vec3 ro = false ? vec3(cos(iTime),0.5,sin(iTime))*10. : texelFetch(BUF_A,RO_COO,0).xyz;
        vec3 lookat = vec3(0);
        mat3 cam = camera(ro, lookat, 0.);
        vec3 rd = cam * normalize(vec3(uv,1));
        float rad = STAR_RAD;
        vec2 sph = asphere(ro,rd,rad);
        float d=  sph.x, d_refl = sph.y;
        vec3 ero = ro+rd*d+d*0.03, enorm = normalize(ero), erd = reflect(rd,enorm);
        erd = refract(erd,enorm,1.);
        HitInfo hit = intersect(ero,erd);
        vec3 col = vec3(0);
        if(hit.dist < 40. && sph.x < 20.)
        {
            //hit.nor = normal(ro + rd*hit.dist);
            hit.nor = normal(ero+erd*hit.dist);
            vec3 lig_pos = vec3(1,10,0);
            col = calcLight(hit,erd,lig_pos);    
            //col *= texelFetch(BUF_A,ivec2(hit.id,1),0).xxx*2.+0.1; 
            col = max(col,vec3(0));        
        }
        vec3 l = normalize(vec3(.7, 0.8, 0));
        vec3 fCol = getSky(rd, l);
        //if(sph.x < MAX_DIST) col = mix(clamp(col, 0., 1.), fCol, smoothstep(.14, .6, sph.y/MAX_DIST));
        //if(false) 
        {
            
            float d_refl = asphere(ro,rd,rad).y; //intersectSphere(ro,rd,vec3(0),rad);
            if(d_refl < 100.)
            {
                vec3 rro = ro+rd*d_refl+d_refl*0.03, rnorm = normalize(rro), rrd = reflect(rd,rnorm);
                rrd = refract(rrd,rnorm,.55);
                float sph_diff = max(0,dot(-rnorm,vec3(2,5,0)))*0.5+0.5;
                vec3 c_sph = vec3(0.9059, 0.9529, 0.0627)*sph_diff*0.01;
                float fres = fresnel(.2,1.14,3.,rro-ro,-rnorm);
                //col += max(vec3(0),c_sph+vec3(fres*0.02));
                HitInfo reflection = intersect(rro,rrd);
                if(reflection.dist < MAX_DIST)
                {
                    //hit.nor = normal(ro + rd*hit.dist);
                    reflection.nor = normal(rro+rrd*reflection.dist);
                    vec3 lig_pos = vec3(1,0,0);
                    vec3 rcol = calcLight(reflection,rrd,lig_pos);
                    //col *= texelFetch(BUF_A,ivec2(hit.id,1),0).xxx*2.+0.1; 
                    rcol = max(rcol,vec3(0));
                    col = mix(rcol,col,0.5);        
                }

            }
           
        
        }
        
        vec3 feed = max(pow(texture(FEEDBACK,fragCoord/iResolution.xy).xyz,vec3(2.)),vec3(0));
        //col = mix(col, feed, 0.2)*0.051+col*0.25;
        int song_sect = getSongSection(MIDI);
        vec3 sun1_pos = perhi_sect_displ(song_sect);
        vec3 sun1 = is_perhi_on(song_sect) ? integrateLightFullView(ro-sun1_pos,rd,0.051,1.)*PERHI_COL_CENTER : vec3(0);
        vec3 sun2_pos = flower_sect_displ(song_sect);
        vec3 sun2 = is_flower_on(song_sect) ?  integrateLightFullView(ro-sun2_pos,rd,0.21,0.5)*FLOWER_COL_CENTER :  vec3(0);
        //col += sun1+sun2*0.1;
        
        col = pow(col,vec3(.6445));
        tot += col;
#if AA > 1
    }

    tot /= float(AA);
#endif
    fragColor= = vec4(pow(tot, vec3(.4545),glow));
}


