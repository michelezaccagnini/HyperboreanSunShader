#define FEEDBACK iChannel0
#define BUF_A iChannel1
#define MIDI iChannel2
#define TEXTURE iChannel3

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
        res = map(p, BUF_A, MIDI, TEXTURE);
        if(d > max_dist || abs(res.dist) < min_dist) break;
        d += res.dist;
    }
    res.dist = d;
    return res;
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
    stripe = pow(smoothstep(0.1,0.0,abs(pow(hit.env,0.5)-(1.-hit.uv.y))),2.);
    col *=vec3(stripe *3.+10.8);//diff*pow(hit.env,1.5)*1.+0.1;//max(col, vec3(0));
    vec3 txt = texture(TEXTURE,hit.uv_transorm).xyz;
    float bump = clamp(dot(txt,txt),0.,1.);
    //col += txt*light_intensity;
    col += bump*PERHI_COL_CENTER;
    // color center spheres for FLOWER and PERHI
    col =  mix(col,sun_col,hit.smin);
    if(hit.id.x == 1 && hit.id.y > 4) col *= 0.2;
    return col;
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
    vec3 col = vec3(0);
    if(hit.dist < 40.)
    {
        //hit.nor = normal(ro + rd*hit.dist);
        hit.nor = normal(ro+rd*hit.dist, BUF_A, MIDI, TEXTURE);
        vec3 lig_pos = vec3(1,10,0);
        col = calcLight(hit,rd,lig_pos);
        //col *= texelFetch(BUF_A,ivec2(hit.id,1),0).xxx*2.+0.1; 
        col = max(col,vec3(0));        
    }
    vec3 feed = max(pow(texture(FEEDBACK,fragCoord/iResolution.xy).xyz,vec3(2.)),vec3(0));
    col = mix(col, feed, 0.2)*0.5+col*0.8;
    int song_sect = getSongSection(MIDI);
    vec3 sun1_pos = perhi_sect_displ(song_sect);
    vec3 sun1 = is_perhi_on(song_sect) ? integrateLightFullView(ro-sun1_pos,rd,0.051,1.)*PERHI_COL_CENTER : vec3(0);
    vec3 sun2_pos = flower_sect_displ(song_sect);
    vec3 sun2 = is_flower_on(song_sect) ?  integrateLightFullView(ro-sun2_pos,rd,0.21,0.5)*FLOWER_COL_CENTER :  vec3(0);
    fragColor.xyz = pow(col+sun1+sun2, vec3(.4545));
    // fragColor = texelFetch(iChannel2, ivec2(fragCoord.xy/iResolution.xy*vec2(128.,16.*5.)), 0);
}


