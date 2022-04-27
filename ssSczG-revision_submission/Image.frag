#define BUF_B iChannel1
#define BUF_A iChannel0
#define BUF_C iChannel2
#define MIDI iChannel3
#define BUMPFACTOR 0.01

vec3 get_uv(vec3 p)
{
    float rep =2., sca = 0.8;
    vec2 uv = polar(p)*rep;
    uv.x = mod(uv.x,1.)*0.65+0.25, uv.y = mod(uv.y,1.)*0.65+0.25;
    float vignette = smoothstep(0.9,0.6,length(uv));
    return vec3(uv,vignette);
}

float get_bump(vec3 p)
{
    vec3 uv = get_uv(p);
    vec3 t = texture(BUF_B,uv.xy).xyz*0.8*uv.z;
    return t.x*0.5;
}
vec3 normal(vec3 ro, vec3 rd, float rad, inout vec3 norm_back)
{
    float dx = 5.e-4;
    vec2 c = vec2(1,0);
    vec2 s = asphere(ro,rd,rad);
   
    float bumf= get_bump(ro+rd*s.x)*BUMPFACTOR, bumb = get_bump(ro+rd*s.y)*BUMPFACTOR;
    norm_back =normalize(vec3(
                          asphere(ro,rd+dx*c.xyy,rad-bumb).y,
                          asphere(ro,rd+dx*c.yxy,rad-bumb).y,
                          asphere(ro,rd+dx*c.yyx,rad-bumb).y
                        ) - s.y);
    return normalize(vec3(
              asphere(ro,rd+dx*c.xyy,rad-bumf).x,
              asphere(ro,rd+dx*c.yxy,rad-bumf).x,
              asphere(ro,rd+dx*c.yyx,rad-bumf).x
            ) - s.x);
}
vec3 blinn_phong(vec3 p, vec3 rd, vec3 light, vec3 norm)
{
    vec3 col_ambient = vec3(0.1);
    vec3 col_diffuse =vec3(0.5,0.5,0.1);
    vec3 col_specular = vec3(0.1,0.2,0.85);
    return  col_ambient + 
            col_diffuse * max(dot(normalize(light-p),norm),0.)+ 
            col_specular * pow(max(dot(reflect(normalize(light-p),norm),rd),0.),2.);

}

float sphDensity( vec3  ro, vec3  rd,   // ray origin, ray direction
                  vec3  sc, float sr,   // sphere center, sphere radius
                  float dbuffer )       // depth buffer
{
    // normalize the problem to the canonical sphere
    float ndbuffer = dbuffer / sr;
    vec3  rc = (ro - sc)/sr;
	
    // find intersection with sphere
    float b = dot(rd,rc);
    float c = dot(rc,rc) - 1.0;
    float h = b*b - c;

    // not intersecting
    if( h<0.0 ) return 0.0;
	
    h = sqrt( h );
    
    //return h*h*h;

    float t1 = -b - h;
    float t2 = -b + h;

    // not visible (behind camera or behind ndbuffer)
    if( t2<0.0 || t1>ndbuffer ) return 0.0;

    // clip integration segment from camera to ndbuffer
    t1 = max( t1, 0.0 );
    t2 = min( t2, ndbuffer );

    // analytical integration of an inverse squared density
    float i1 = -(c*t1 + b*t1*t1 + t1*t1*t1/3.0);
    float i2 = -(c*t2 + b*t2*t2 + t2*t2*t2/3.0);
    return (i2-i1)*(3.0/4.0);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = (2.0*(fragCoord)-iResolution.xy)/iResolution.y;
    
    vec4 bufb = texture(BUF_B,fragCoord/iResolution.xy);//max(texture(BUF_C,fragCoord/iResolution.xy).xyz,texture(BUF_B,fragCoord/iResolution.xy).xyz); 
    vec3 col = bufb.xyz, glow = FLOWER_COL_CENTER*bufb.w;
    vec3 ro = false ? vec3(cos(iTime),0.5,sin(iTime))*10. : texelFetch(BUF_A,RO_COO,0).xyz;
    vec3 lookat = vec3(0);
    mat3 cam = camera(ro, lookat, 0.);
    vec3 rd = cam * normalize(vec3(uv,1));
    float rad = STAR_RAD;
    vec2 sph = asphere(ro,rd,STAR_RAD);
    float b = get_bump(ro+rd*sph.x);
    //sph.x = asphere(ro,rd,rad-b).x;
    //b = get_bump(ro+rd*sph.y);
    //sph.y = asphere(ro,rd,rad-b).y;
    float front=  sph.x, back = sph.y;
    vec3 fr = get_uv(ro+rd*front), bk = get_uv(ro+rd*back);
    vec3 txf = texture(BUF_B,fr.xy).xyz*1.0*fr.z;
    vec3 txb = texture(BUF_B,bk.xy).xyz*1.0*bk.z;
    float dens = sphDensity(ro,rd,vec3(0),STAR_RAD,1000.);
    vec3 p = ro+rd*sph.x;
    
    if(sph.x < MAX_DIST && sph.x > 0.2)
    {
        glow *= 0.1;
        vec3 nba = vec3(0);
        vec3 nfr = normal(ro,rd,STAR_RAD,nba);
        col += blinn_phong(ro+rd*front, rd, vec3(0,1,0),nfr)*(glow+0.1);
        col += blinn_phong(ro+rd*back, rd, vec3(0,1,0),nba)*(glow+0.1);
        col += txf+txb;
        //col = nfr;
        //col = vec3(sph_uv1.xxx);
    }


    /*
    vec2 o = vec2(0.0015,0);
    col += texture(BUF_B,fragCoord/iResolution.xy+o).xyz;
    col += texture(BUF_B,fragCoord/iResolution.xy+o.yx).xyz;
    col += texture(BUF_B,fragCoord/iResolution.xy-o).xyz;
    col += texture(BUF_B,fragCoord/iResolution.xy-o.yx).xyz;
    col /= 4.;
    */
    #if 0
    //fragColor = vec4(getSongSection(MIDI) == 2 ? 1 : 0);
    
    fragColor = false ? texelFetch(BUF_C,ivec2(fragCoord*vec2(1)),0) : 
                    texelFetch(BUF_A,ivec2(fragCoord*vec2(0.03,0.251)),0).xxxx  ;
    
    #else
    //col = dot(col,col) > 0.01 ? min(vec3(glow),col): col;
    //glow *= dot(col,col) > 0.1 ? 0. : 2.;
    fragColor = vec4(pow(col+vec3(clamp(glow,0.,1.)),vec3(.6545)),0.);
    #endif
}