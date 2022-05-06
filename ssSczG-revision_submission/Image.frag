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




void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = (2.0*(fragCoord)-iResolution.xy)/iResolution.y;
    
    vec4 bufb = texture(BUF_B,fragCoord/iResolution.xy);//max(texture(BUF_C,fragCoord/iResolution.xy).xyz,texture(BUF_B,fragCoord/iResolution.xy).xyz); 
    vec3 col = bufb.xyz, glow = FLOWER_COL_CENTER*bufb.w;
    float fgl = bufb.w ;
    vec3 ro = false ? vec3(cos(iTime),0.5,sin(iTime))*10. : texelFetch(BUF_A,RO_COO,0).xyz;
    vec3 lookat = LOOKAT;
    mat3 cam = camera(ro, lookat, 0.);
    vec3 rd = cam * normalize(vec3(uv,1));
    float rad = STAR_RAD;
    vec2 sph = asphere(ro,rd,STAR_RAD);
    float b = get_bump(ro+rd*sph.x);
    //sph.x = asphere(ro,rd,rad-b).x;
    //b = get_bump(ro+rd*sph.y);
    //sph.y = asphere(ro,rd,rad-b).y;
    float front=  sph.x, back = sph.y;
    vec3 uv_fr = get_uv(ro+rd*front), uv_bk = get_uv(ro+rd*back);
    vec3 txf = texture(ORGA,uv_fr.xy).xyz*1.0*uv_fr.z;
    vec3 txb = texture(ORGA,uv_bk.xy).xyz*1.0*uv_bk.z;
    vec4 dru = texture(BUF_C,fragCoord/iResolution.xy);
    fgl -= dru.w;
    if(sph.x < MAX_DIST && sph.x > 0.2)
    {
        //glow*= 0.1;
        fgl = smoothstep(2.2,3.9,fgl*30.)*0.5;//pow(clamp(fgl*0.5,0.,1.),2.)*0.5;
        vec3 norm_front = normalize(ro+rd*sph.x), norm_back = normalize(ro+rd*sph.y);
        // bump mapping
	    vec3 surf_norm_front = get_bump_norm(norm_front,uv_fr.xy), 
             surf_norm_back = get_bump_norm(norm_back , uv_bk.xy);
        
        col = max(blinn_phong(ro+rd*front, rd, vec3(0,1,0),surf_norm_front, SunCol)*fgl,col);
        col = max(blinn_phong(ro+rd*back , rd, vec3(0,1,0),surf_norm_back, SunCol )*fgl,col);
        col += (txf+txb*0.5)*(fgl+dru.w)*SunCol;
        //col = vec3(glow);
        //col = nfr;
        //col = vec3(sph_uv1.xxx);
    }


    #if 0
    vec2 o = vec2(0.0015,0);
    col += texture(BUF_B,fragCoord/iResolution.xy+o).xyz;
    col += texture(BUF_B,fragCoord/iResolution.xy+o.yx).xyz;
    col += texture(BUF_B,fragCoord/iResolution.xy-o).xyz;
    col += texture(BUF_B,fragCoord/iResolution.xy-o.yx).xyz;
    col /= 4.;
    #endif
    col = pow(col+clamp(glow,0.,1.),vec3(.4545));
    
    #if 0
    //fragColor = vec4(getSongSection(MIDI) == 2 ? 1 : 0);
    
    fragColor = true ?  texelFetch(BUF_C,ivec2(fragCoord),0) : vec4(1);
    
    #else
    
    //col = dot(col,col) > 0.01 ? min(vec3(glow),col): col;
    //glow *= dot(col,col) > 0.1 ? 0. : 2.;
    fragColor = vec4(col+dru.xyz,0.);
    #endif
}