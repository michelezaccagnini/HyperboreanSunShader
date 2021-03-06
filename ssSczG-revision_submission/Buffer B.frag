//This is the main buffer where all the objects are mapped
// the glow value (moving stars and local glow) is passed on the alpha chan

#define FEEDBACK iChannel0
#define BUF_A iChannel1
#define MIDI iChannel2
#define TEXTURE iChannel3
#define ZERO (min(iFrame,0))

float glow;
float glow_perhi;
float glow_flower;
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
    vec3 hit_point = vec3(0);
    //moving star dust
    vec3 op = p;
    op.xy *= rotate(iTime*0.12);
    op.z += iTime*ASTEROID_SPEED;
    float pa =  length( cos(op*.7)+cos(op*0.812)*0.1);
    glow +=.18/(0.1+pa*pa*pa*500.);
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
            float c_sphere = length(p-sp_pos)-0.5;
            hit_point = hit_point;
            float h = 0.;
            float d = rop.x;
            rop.x = smin(c_sphere,d,0.3,h);
            float o_sphere = length(p)-STAR_RAD*0.9;
            rop.x = smax(rop.x,o_sphere,0.3);
            float glow_int = smoothstep(0.051,0.,abs(rop.z*1.-(1.-pow(env,.75))))*smoothstep(STAR_RAD,STAR_RAD/2.,length(p));
            glow += (.1/(0.1+rop.x*rop.x*rop.x))*glow_int;   
            glow_flower += (0.01/(0.1+rop.x*rop.x*rop.x*rop.x+20.));
            
            if(c_sphere< d && !refl)
            {
                glow += smoothstep(0.3,0.2,(0.01/(0.1+rop.x*rop.x+cos(p.z)*80.))*pow(h,0.5))*0.0015;     
            } 
            if(rop.x < res.dist)
            {
                //glow = min(glow_flower,glow);
                //id += is_first ? 0 : 5;
                vec2 tuv = rop.yz*0.3;
                tuv.y += env;
                vec3 txt = texture(TEXTURE,tuv).xyz;
                float bump = clamp(txt.r*0.1,0.,1.)*0.41;
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
            float glow_int = smoothstep(0.051,0.,abs(rop.z*0.5-env))*1.5*env;
            glow += (.2/(0.1+rop.x*rop.x*rop.x))*glow_int;
            //rop.x = min(rop.x,sph);
            if(rop.x < res.dist)
            {
                vec2 tuv = rop.yz*0.1;
                tuv.y += env;
                vec3 txt = texture(TEXTURE,tuv).xyz;
                float bump = clamp(txt.r,0.,1.)*.21;
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
            float c_sphere = length(p-sp_pos)-.2;
            hit_point = hit_point;
            float h = 0.;
            float d = rop.x;
            rop.x = smin(c_sphere,d,0.9,h);
            float o_sphere = length(p)-STAR_RAD*0.9;
            rop.x = smax(rop.x,o_sphere,0.5);
            float glow_int = smoothstep(0.025,0.02,abs(rop.z*0.2-pow(env,0.6)))*pow(rop.z,1.5)*1.;//*smoothstep(0.5,0.4,pow(rop.z,2.5));
            glow += (0.1/(0.1+rop.x*rop.x))*glow_int;
            glow_perhi += (0.01/(0.1+rop.x*rop.x*rop.x*rop.x+20.));
            if(c_sphere< d && !refl)
            {
                //glow += smoothstep(0.3,0.2,(0.1/(0.1+rop.x*rop.x))*pow(1.-h,0.5))*0.0015;
                
            } 
            if(rop.x < res.dist)
            {
                //glow = min(glow_perhi,glow);
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
   //drums and bass are rendered in buf C and D as lights
    
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
    float glo_sli = .6;
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
            if(hit.id.x == 2)
            {
                col = get_lightPong(hit, rd, l, FLOWER_COL_CENTER, TEXTURE);
            } 
            else if(hit.id.x == 3)
            {
                col = get_lightPerhi(hit, rd, l, FLOWER_COL_CENTER, TEXTURE);
                col = mix(col,SunCol*0.2, 1.-hit.smin);
            } 

            else
            {
                col = get_lightFlower(hit, rd, l, FLOWER_COL_CENTER, TEXTURE);
                col = mix(col,SunCol*0.2, hit.smin);
            }   
            glow = 0.;  
            occl = 0.;   
        }
        vec3 pla = plane(ro,rd,vec3(0),normalize(ro));
        float pl_dist = max(0.,mapSimple(pla));
        float d_ctr = sphDensity(ro,rd,vec3(0),STAR_RAD*2,10000.);
        vec3 feed = max(pow(texture(FEEDBACK,fragCoord/iResolution.xy).xyz,vec3(2.)),vec3(0));
        int song_sect = getSongSection(MIDI);
        vec3 sun =  clamp(integrateLightFullView(ro-SunPos,rd,5.81*pl_dist,0.18),0.,1.)*SunCol;
        col += pow(sun,vec3(2.5))*occl*0.4/float(AA);

        
        tot += max(col,vec3(0));

#if AA > 1
    }
    tot /= float(AA);
    glow /= float(AA*3);
    glow_perhi /= float(AA);
    glow_flower /= float(AA);
#endif
    tot = pow(tot,vec3(.4545));
    float g = slide(texture(FEEDBACK,fragCoord/iResolution.xy).w,glow,glo_sli*STREAM_SLIDE)*.4+pow(glow,1.5)*0.0051;
    g = max(mix(glow_perhi,glow_flower,0.5)*1.,g);
    fragColor= vec4(tot, g);
}


