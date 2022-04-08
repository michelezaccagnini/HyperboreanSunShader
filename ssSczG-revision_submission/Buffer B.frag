#define FEEDBACK iChannel0
#define BUF_A iChannel1
#define MIDI iChannel2
#define TEXTURE iChannel3


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
    vec3 pos;
    vec3 surf;
    vec3 nor;
    vec2 uv;
    vec3 col;
    float env;
};

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
    vec3 hit_point = vec3(0), norm = vec3(0);
    //flower
    if(flower_on)
    {
        vec3 hit_point = vec3(0),
             hit_point2 = vec3(0),
             norm = vec3(0), 
             norm2 = vec3(0);
        for(int id = 0; id < 5; id++)
        {
            
            vec3 rope_center = rope_flower1(p,id+FLOWER_BLOCK.z  ,BUF_A,0.01,hit_point,norm);
            vec3 rope_attach = rope_flower2(p,id+FLOWER_BLOCK.z*3,BUF_A,0.01,hit_point2,norm2); 
            //float joint = length(p-texelFetch(BUF_A,ivec2(0,id+FLOWER_BLOCK.z*3),0).xyz)-0.1;
            vec3 rop = opU(rope_center,rope_attach);
            bool is_first = rope_center.x < rope_attach.x; 
            hit_point = is_first ? hit_point : hit_point2;
            norm = is_first ? norm : norm2 ;
            //rop.x = smin(joint,rop.x,0.7);
            if(rop.x < res.dist)
            {
                res.dist = rop.x, res.id = ivec2(1,id), res.uv = rop.yz, 
                res.pos = p, res.surf = hit_point,  res.nor = norm;
                res.env = texelFetch(BUF_A,ivec2(id,FLOWER_ENV_ROW),0).x;            
            }
        }
    }
    if(pong_on)
    {
        vec3 hit_point = vec3(0), norm = vec3(0);
        for(int id = 0; id < 8; id++)
        {
            float env = texelFetch(BUF_A,ivec2(id,PONG_ENV_ROW),0).x;
            vec3 pos = texelFetch(BUF_A,ivec2(0,id+PONG_BLOCK_OFFSET+PONG_BLOCK.z),0).xyz;
            //float sph = length(p - pos) - 0.4;
            vec3 rop  = rope(p, id+PONG_BLOCK_OFFSET+PONG_BLOCK.z,BUF_A,env, hit_point, norm);
            //rop.x = min(rop.x,sph);
            if(rop.x < res.dist)
            {
                res.dist = rop.x, res.id = ivec2(2,id), 
                res.uv = rop.yz, res.pos = p, res.surf = hit_point, 
                res.env = texelFetch(BUF_A,ivec2(id,PONG_ENV_ROW),0).x;
                res.surf = hit_point, res.nor = norm;
            }
        }
    }
    if(perhi_on )
    {
        vec3 hit_point = vec3(0), norm = vec3(0);
        for(int id = 0; id < 8; id++)
        {
            float env = texelFetch(BUF_A,ivec2(id,PERHI_ENV_ROW),0).x;
            if(env< 0.01) continue;
            vec3 pos = texelFetch(BUF_A,ivec2(0,id+PERHI_BLOCK_OFFSET+PERHI_BLOCK.z),0).xyz;
            float sph = length(p - pos) - 0.4;
            float center_sph = length(p) -0.4;
            vec3 rop  = ropePerhi(p, id+PERHI_BLOCK_OFFSET+PERHI_BLOCK.z,BUF_A,0.01, env, hit_point, norm);
            //rop.x = min(rop.x,sph);
            //rop.x = min(rop.x,center_sph);
            if(rop.x < res.dist)
            {
                res.dist = rop.x, res.id = ivec2(3,id), res.uv = rop.yz, res.pos = p, 
                res.surf = hit_point, 
                res.env = texelFetch(BUF_A,ivec2(id,PERHI_ENV_ROW),0).x;
                res.surf = hit_point, res.nor = norm;
            }
        }
    }
    if(drums_on)
    {
        vec3 hit_point = vec3(0), norm = vec3(0);
        for(int id = 0; id < 8; id++)
        {
            vec3 pos = texelFetch(BUF_A,ivec2(0,id+DRUMS_BLOCK_OFFSET),0).xyz;
            float env = texelFetch(BUF_A,ivec2(id,DRUMS_ENV_ROW),0).x;
            float sph = length(p - pos) - 0.2;
            vec3 rop  = ropeDrums(p, id+DRUMS_BLOCK_OFFSET+DRUMS_BLOCK.z,BUF_A,0.01, env, hit_point, norm);
            //rop.x = min(rop.x,sph);
            if(rop.x < res.dist)
            {
                res.dist = rop.x, res.id = ivec2(4,id), res.uv = rop.yz, 
                res.pos = p, res.surf = hit_point, 
                res.env = texelFetch(BUF_A,ivec2(id,DRUMS_ENV_ROW),0).x;
                res.surf = hit_point, res.nor = norm;

            }
        }
    }
    if(bass_on)
    {
        vec3 hit_point = vec3(0), norm = vec3(0);
        for(int id = 0; id < 8; id++)
        {
            vec3 pos = texelFetch(BUF_A,ivec2(0,id+BASS_BLOCK_OFFSET),0).xyz;
            float sph = length(p - pos) - 0.2;
            vec3 rop  = rope(p, id+BASS_BLOCK_OFFSET+BASS_BLOCK.z,BUF_A,0.01, hit_point, norm);
            //rop.x = min(rop.x,sph);
            if(rop.x < res.dist)
            {
                res.dist = rop.x, res.id = ivec2(4,id), res.uv = rop.yz, 
                res.pos = p, res.surf = hit_point, 
                res.env = texelFetch(BUF_A,ivec2(id,BASS_ENV_ROW),0).x;
                res.surf = hit_point, res.nor = norm;
            }
        }
    }
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
    return res;
}


// http://iquilezles.org/www/articles/normalsSDF/normalsSDF.htm
vec3 normal( in vec3 pos)
{
    const float ep = 0.01;
    vec2 e = vec2(1.0,-1.0)*0.05773;
    return normalize( e.xyy*map( pos + e.xyy*ep).dist + 
					  e.yyx*map( pos + e.yyx*ep).dist + 
					  e.yxy*map( pos + e.yxy*ep).dist + 
					  e.xxx*map( pos + e.xxx*ep).dist );
}

vec3 calcLight(HitInfo hit, vec3 rd, vec3 lig_pos)
{ 
    vec3 diff = max(0.,dot(hit.surf, lig_pos)*0.5+0.5)*vec3(0.2,1.12,1.26);
    //reflect
    vec3 refl = reflect(lig_pos,hit.pos);
    vec3 spec_lig = vec3(clamp(dot(rd, refl),0.,1.));
    float light_distance = smoothstep(.1,1.01, length(lig_pos- hit.pos)*.1)*1.2;
    float falloff = .15;
    float light_intensity = clamp(falloff/pow(light_distance,0.42),0,1);
    //light_intensity =pow(light_intensity,.4545);
    float fres = clamp(1. - dot(hit.nor, -rd),0.,1.);
    vec3 col = fres*0.2*diff*vec3(pow(hit.env,.5)*5.)*hit.col*1.;
    col *=vec3(hit.env *3.+0.1);//diff*pow(hit.env,1.5)*1.+0.1;//max(col, vec3(0));
    vec3 txt = texture(TEXTURE,hit.uv*0.01).xyz;
    return mix(col,max(col,txt),0.4);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Normalized pixel coordinates (from 0 to 1)
    vec2 uv = (fragCoord-0.5*iResolution.xy)/iResolution.y*2.;
    vec3 ro = texelFetch(BUF_A,RO_COO,0).xyz;
    //ro.xz *= rotate(mod(iTime,3.36)/3.36*TAU);
    vec3 lookat = vec3(0);
    mat3 cam = camera(ro, lookat, 0.);
    vec3 rd = cam * normalize(vec3(uv,1));
    HitInfo hit = intersect(ro,rd);
    hit.col =  hash31(float(hit.id.x*8)*float(hit.id.y)*0.3);
    //hit.env = texelFetch(BUF_A,ivec2(hit.id,FLOWER_ENV_ROW),0).x;
    vec3 col = vec3(0);
    if(hit.dist < 10.)
    {
        //hit.nor = normal(ro + rd*hit.dist);
        vec3 lig_pos = vec3(1,10,0);
        col = calcLight(hit,rd,lig_pos);
        //col *= texelFetch(BUF_A,ivec2(hit.id,1),0).xxx*2.+0.1; 
        col = max(col,vec3(0));
        
        //col = mix(col, feed,0.1);
    }
    //add PERHI textured sphere
    //if(sph.x < 100.) col += (txt1+txt2);
    //col =txt1 + txt2;
    //fragColor = texelFetch(BUF_A,ivec2(0.02*fragCoord),0);
    fragColor.xyz = pow(col, vec3(.4545));
}


