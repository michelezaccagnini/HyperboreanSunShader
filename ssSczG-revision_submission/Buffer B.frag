#define FEEDBACK iChannel0
#define BUF_A iChannel1
#define MIDI iChannel2
#define TEXTURE iChannel3
#define ZERO (min(iFrame,0))




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

float glow;
#define SUN_COLOR vec3(0.8902, 0.9333, 0.0784)
#define SUN_POS vec3(0,0,100)
vec3 sun;

float mapSimple(vec3 p)
{
    float dist = MAX_DIST;
    int song_sect = getSongSection(MIDI);
    bool flower_on = is_flower_on(song_sect);
    bool pong_on = is_pong_on(song_sect);
    bool perhi_on = is_perhi_on(song_sect);
    bool drums_on = is_drums_on(song_sect);
    bool bass_on = is_bass_on(song_sect);
    vec3 hit_point;
    if(flower_on)
    {
        for(int id = 0; id < 5; id++)
        {
            vec3 displ = flower_sect_displ(song_sect);
            vec3 sp_pos = displ;
            float env = texelFetch(BUF_A,ivec2(id,FLOWER_ENV_ROW),0).x;
            vec3 rop = rope_flower1(p-displ,id+FLOWER_BLOCK.z  ,BUF_A,env,hit_point);
            float c_sphere = length(p-sp_pos)-0.5;
            float h = 0.;
            float d = rop.x;
            rop.x = smin(c_sphere,d,0.3,h);
            if(rop.x < dist) dist = rop.x;
        }
    }
    if(pong_on)
    {
        for(int id = 0; id < 8; id++)
        {
            float env = texelFetch(BUF_A,ivec2(id,PONG_ENV_ROW),0).x;
            vec3 rop  = ropePong(p, id+PONG_BLOCK_OFFSET+PONG_BLOCK.z,BUF_A,env, hit_point);
            if(rop.x < dist)  dist = rop.x;
        }
    }
    if(perhi_on )
    {
        ivec2 o_sect = any(equal(ivec2(song_sect),ivec2(2,4))) ? ivec2(0 ,4) : ivec2(0,8) ; 
        for(int id = o_sect.x; id < o_sect.y; id++)
        {
            float env = texelFetch(BUF_A,ivec2(id,PERHI_ENV_ROW),0).x;
            vec3 displ = perhi_sect_displ(song_sect);
            vec3 sp_pos = displ ;
            vec3 rop  = ropePerhi(p-displ, id+PERHI_BLOCK_OFFSET+PERHI_BLOCK.z,BUF_A, env, hit_point);
            float c_sphere = length(p-sp_pos)-.7;
            float h = 0.;
            float d = rop.x;
            rop.x = smin(c_sphere,d,0.9,h);
            float o_sphere = length(p)-STAR_RAD*2.9;
            rop.x = smax(rop.x,o_sphere,0.5);
            if(rop.x < dist) dist = rop.x; 
        }
    }
    if(drums_on)
    {
        int block = DRUMS_BLOCK_OFFSET+DRUMS_BLOCK.z;
        for(int id = 0; id < 8; id++)
        {
            float env = texelFetch(BUF_A,ivec2(id,DRUMS_ENV_ROW),0).x;
            vec3 rop  = ropeDrums(p, id+block,BUF_A, env, hit_point);
            if(rop.x < dist) dist = rop.x;
        }
    }
    if(bass_on)
    {
        for(int id = 0; id < 8; id++)
        {
            float env = texelFetch(BUF_A,ivec2(id,BASS_ENV_ROW),0).x;
            vec3 displ = song_sect > 5 ?vec3(0,0, LAST_SECT_DISPLACE) : vec3(0);
            vec3 rop  = ropeDrums(p-displ, id+BASS_BLOCK_OFFSET+BASS_BLOCK.z,BUF_A, env, hit_point);
            if(rop.x < dist) dist = rop.x;
        }
    }
    return dist;
}



HitInfo map(vec3 p, bool refl)
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
    vec3 op = p;
    op.z += iTime*10.;
    float b = cos(op.z*.2); ;
    float pa =  length( cos(op*1.7));//max(length(cos(op*.75+vec3(b,b*.5,cos(p.x)))),abs(p.x)-10.);
    glow +=1.8/(0.1+pa*pa*500.);
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
            float c_sphere = length(p-sp_pos)-0.5;
            hit_point = hit_point;
            float h = 0.;
            float d = rop.x;
            rop.x = smin(c_sphere,d,0.3,h);
            float o_sphere = length(p)-STAR_RAD*0.9;
            rop.x = smax(rop.x,o_sphere,0.3);
            float glow_int = smoothstep(0.051,0.,abs(rop.z*1.-(1.-pow(env,.75))))*smoothstep(STAR_RAD,STAR_RAD/2.,length(p));
            glow += (.1/(0.1+rop.x*rop.x*rop.x))*glow_int;
            
            if(c_sphere< d && !refl)
            {
                glow += smoothstep(0.3,0.2,(0.01/(0.1+rop.x*rop.x+cos(p.z)*80.))*pow(h,0.5))*0.0005;
                
            } 
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
            glow += (2.1/(0.1+rop.x*rop.x))*glow_int;
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
            vec3 rop  = ropePerhi(p-displ, id+PERHI_BLOCK_OFFSET+PERHI_BLOCK.z,BUF_A, env, hit_point);
            float c_sphere = length(p-sp_pos)-.7;
            hit_point = hit_point;
            float h = 0.;
            float d = rop.x;
            rop.x = smin(c_sphere,d,0.9,h);
            float o_sphere = length(p)-STAR_RAD*2.9;
            rop.x = smax(rop.x,o_sphere,0.5);
            float glow_int = smoothstep(0.025,0.02,abs(rop.z*0.2-pow(env,0.6)))*pow(rop.z,1.5)*1.;//*smoothstep(0.5,0.4,pow(rop.z,2.5));
            glow += (0.1/(0.1+rop.x*rop.x))*glow_int;
            if(c_sphere< d && !refl)
            {
                glow += smoothstep(0.3,0.2,(0.01/(0.1+rop.x*rop.x))*pow(1.-h,0.5))*0.0005;
                
            } 
            if(rop.x < res.dist)
            {
                vec2 tuv = rop.yz*0.14;
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
            //float d = rop.x;
            //float h = 0.;
            //vec3 displ = perhi_sect_displ(song_sect);
            //vec3 sp_pos = displ ;
            //float c_sphere = length(p-sp_pos)-STAR_RAD*0.8;
            //rop.x = smin(c_sphere,d,0.5,h);
            //float o_sphere = length(p)-STAR_RAD*2.9;
            //rop.x = smax(rop.x,o_sphere,0.5);
            if(rop.x < res.dist)// && h < 0.99999)
            {
                vec2 tuv = rop.yz*0.1;
                env = pow(env,0.5);
                tuv.y += env;
                vec3 txt = texture(TEXTURE,tuv).xyz;
                float bump = clamp(dot(txt,txt),0.,1.)*0.4*(env+0.05);
                float glow_int = smoothstep(0.04,0.,abs(rop.z*0.8-pow(env,2.)))*.5*env;
                glow += (0.1/(0.1+rop.x*rop.x*rop.x))*glow_int;
                res.dist = rop.x-bump, res.id = ivec2(4,id), res.uv_transorm = tuv,
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
    
    return res;
}

HitInfo intersect(vec3 ro, vec3 rd, bool refl)
{
    HitInfo res;
    float min_dist = 0.01;
    float max_dist = MAX_DIST;
    float d = min_dist;
    res.id = ivec2(-1);
    for(int i = 0; i < 80; i++)
    {
        vec3 hit_point = vec3(0);
        
        vec3 p = ro + rd*d;
        res = map(p, refl);
        if(d > max_dist || abs(res.dist) < min_dist) break;
        d += res.dist;
    }
    res.dist = d;
    return res;
}


vec3 normal( in vec3 pos)
{
    const float ep = 5.e-4;
    vec2 e = vec2(1.,0.);
    bool refl = false;
    float s = map(pos,refl).dist;
    
    return normalize( vec3(map( pos + e.xyy*ep,refl).dist, 
                           map( pos + e.yxy*ep,refl).dist, 
                           map( pos + e.yyx*ep,refl).dist) - s);
}



void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec3 tot = vec3(0);
    float glo_sli = 0.5;
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
        vec3 lookat = LOOKAT;
        mat3 cam = camera(ro, lookat, 0.);
        vec3 rd = cam * normalize(vec3(uv,1));
        bool refl = false;
        HitInfo hit = intersect(ro,rd,refl);
        bool is_flo_perhi = any(equal(ivec2(hit.id.x),ivec2(1,3)));
        float rad = STAR_RAD;
        vec2 sph = asphere(ro,rd,rad);
        vec3 l = SunPos;
        vec3 col = vec3(0);
        float occl = 1.;
        if(hit.dist < 40. )
        {   
            glo_sli = 0.4;
            //hit.nor = normal(ro + rd*hit.dist);
            hit.nor = normal(ro+rd*hit.dist);
            //col = calcLight(hit,rd,lig_pos);    
            col = get_light2(hit, rd, l, FLOWER_COL_CENTER, TEXTURE);
            //col *= texelFetch(BUF_A,ivec2(hit.id,1),0).xxx*2.+0.1; 
            col = max(col,vec3(0)); 
            glow = 0.;  
            occl = 0.;   
        }
        vec3 pla = plane(ro,rd,vec3(0),normalize(ro));
        float pl_dist = max(0.,mapSimple(pla));

        
        vec3 fCol = getSky(rd, l);
        //if(sph.x < MAX_DIST) col = mix(clamp(col, 0., 1.), fCol, smoothstep(.14, .6, sph.y/MAX_DIST));
        
        vec3 feed = max(pow(texture(FEEDBACK,fragCoord/iResolution.xy).xyz,vec3(2.)),vec3(0));
        //col = mix(col, feed, 0.2)*0.051+col*0.25;
        int song_sect = getSongSection(MIDI);
        //pl_dist *= smoothstep(8.,3., pl_dist)*smoothstep(100.,8.,pl_dist)+0.2;
        vec3 sun3 =  clamp(integrateLightFullView(ro-SunPos,rd,1.81*pl_dist,0.5),0.,1.)*SunCol;
       
        
        col += pow(sun3,vec3(4.5))*occl;
        
        col = pow(col,vec3(.4545));
        tot += col;
        //tot = vec3(smoothstep(0.,1.,pl_dist*0.1))*0.1;
#if AA > 1
    }

    tot /= float(AA);
    glow /= float(AA);
#endif
    float g = slide(texture(FEEDBACK,fragCoord/iResolution.xy).w,glow,glo_sli*STREAM_SLIDE)*.8*(1.-glo_sli)+pow(glow,2.5)*0.051;
    fragColor= vec4(tot+sun, g);
}


