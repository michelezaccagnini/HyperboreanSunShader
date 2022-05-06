
#define TIME_UNIT 0.12 // 120 ms = 125 bpm
#define FPS 60
const float STREAM_SLIDE = FPS == 60 ? 0.5 : 1.;  
#define MAX_DIST 100.
#define LOOKAT vec3(cos(iTime),sin(iTime),0)*0.8
const vec3 SunPos = vec3(0,0,50);
const vec3 SunCol = vec3(0.9333, 0.8157, 0.1294);

#define AA 1
//#define REDUCE_BOXY

//Flower uniforms
#define ROPE_POINTS 7
#define NUM_FLOWER_PETALS 5
#define FLOWER_BLOCK_OFFSET 0
//block dimensions: x dimension, y dimension, y dim of sub-block, y dimension offset
const ivec4 FLOWER_BLOCK = ivec4(ROPE_POINTS,
                                NUM_FLOWER_PETALS*4+1,//21
                                NUM_FLOWER_PETALS, //5
                                FLOWER_BLOCK_OFFSET);
const int FLOWER_ENV_ROW = FLOWER_BLOCK.y-1;
#define FLOWER_COL_CENTER vec3(0.3059, 0.4039, 0.8941)
const float FLOWER_LENGTH = 4.;

//Pong uniforms
#define PONG_POINTS 7
#define NUM_PONG 8
const int PONG_BLOCK_OFFSET = FLOWER_BLOCK.y;//21
//block dimensions: x dimension, y dimension, y dim of sub-block, y dimension offset
const ivec4 PONG_BLOCK = ivec4(PONG_POINTS,NUM_PONG*2+1, NUM_PONG,PONG_BLOCK_OFFSET); 
const int PONG_ENV_ROW = PONG_BLOCK_OFFSET + PONG_BLOCK.y-1;

//PERHI uniforms
#define PERHI_POINTS 7
#define NUM_PERHI 8
const int PERHI_BLOCK_OFFSET = PONG_BLOCK.y+PONG_BLOCK_OFFSET;//38
//block dimensions: x dimension, y dimension, y dim of sub-block, y dimension offset
const ivec4 PERHI_BLOCK = ivec4(PERHI_POINTS, NUM_PERHI*2+1, NUM_PERHI,PERHI_BLOCK_OFFSET);
const int PERHI_ENV_ROW = PERHI_BLOCK_OFFSET + PERHI_BLOCK.y-1;
#define PERHI_COL_CENTER vec3(0.9451, 0.8549, 0.0392)


#define LAST_SECT_DISPLACE 10.

//PERLO uniforms... PERLO is not used but keep the block uniforms 
//otherwise all other blocks will break :(
#define PERLO_POINTS 7
#define NUM_PERLO 4
const int PERLO_BLOCK_OFFSET = PERHI_BLOCK.y+PERHI_BLOCK_OFFSET;
//block dimensions: x dimension, y dimension, y dim of sub-block, y dimension offset
const ivec4 PERLO_BLOCK = ivec4(PERLO_POINTS, NUM_PERLO*2+1, NUM_PERLO,PERLO_BLOCK_OFFSET);
const int PERLO_ENV_ROW = PERLO_BLOCK_OFFSET + PERLO_BLOCK.y-1;



//DRUMS uniforms
#define DRUMS_POINTS 7
#define NUM_DRUMS 8
#define DRUMS_PITCH_OFFSET 24
#define DRUMS_CHAN 15
const int DRUMS_BLOCK_OFFSET = PERLO_BLOCK.y+PERLO_BLOCK_OFFSET;
//block dimensions: x dimension, y dimension, y dim of sub-block, y dimension offset
const ivec4 DRUMS_BLOCK = ivec4(DRUMS_POINTS, NUM_DRUMS*2+1, NUM_DRUMS,DRUMS_BLOCK_OFFSET);
const int DRUMS_ENV_ROW = DRUMS_BLOCK_OFFSET + DRUMS_BLOCK.y-1;
const int[NUM_DRUMS] DRUMS_TIME_WIDTH = int[NUM_DRUMS](3,4,3,3,4,3,3,3);
#define DRUMS_LENGTH 0.72
#define ASTEROID_SPEED 10.


//Camera uniforms
#define RO_CHAN 14
const int RO_BLOCK_OFFSET = DRUMS_BLOCK.y+DRUMS_BLOCK_OFFSET;
const ivec4 RO_BLOCK = ivec4(1,1,1,RO_BLOCK_OFFSET);
const ivec4 RO_CC = ivec4(80,81,82,83);
const ivec2 RO_COO = ivec2(0,RO_BLOCK_OFFSET);
#define RO_DIST_MULT 50.

//Bass uniforms
#define BASS_POINTS 7
#define NUM_BASS 8
#define BASS_PITCH_OFFSET 24
#define BASS_CHAN 14 
const int BASS_BLOCK_OFFSET = RO_BLOCK.y + RO_BLOCK_OFFSET;
const ivec4 BASS_BLOCK = ivec4(BASS_POINTS, NUM_BASS*2+1, NUM_BASS, BASS_BLOCK_OFFSET);
const int BASS_ENV_ROW = BASS_BLOCK_OFFSET + BASS_BLOCK.y - 1;
#define BASS_TUBE_RADIUS 8.
#define BASS_LENGTH 0.5


//outer sphere
#define STAR_RAD 5.5    

#define TAU 6.28318530718
#define PI 3.14159265359
#define RHO 1.570796326795



mat3 camera( in vec3 ro, in vec3 ta, float cr )
{
	vec3 cw = normalize(ta-ro);
	vec3 cp = vec3(sin(cr), cos(cr),0.0);
	vec3 cu = normalize( cross(cw,cp) );
	vec3 cv = normalize( cross(cu,cw) );
    return mat3( cu, cv, cw );
}

float slide(float cur, float tar, float sli)
{
    return mix(cur, tar, sli);
}
vec2 slide(vec2 cur, vec2 tar, float sli)
{ 
    return mix(cur, tar, sli);
}
vec3 slide(vec3 cur, vec3 tar, float sli)
{ 
    return mix(cur, tar, sli);
}
vec4 slide(vec4 cur, vec4 tar, float sli)
{ 
    return mix(cur, tar, sli);
}

vec3 pix_stream(ivec2 coo, sampler2D text, float sli)
{
    vec3 tar = texelFetch(text, coo - ivec2(1,0),0).xyz;
    vec3 cur = texelFetch(text, coo ,0).xyz;
    return  slide(cur,tar,sli);
}

vec4 pix_stream4(ivec2 coo, sampler2D text, float sli)
{
    vec4 tar = texelFetch(text, coo - ivec2(1,0),0);
    vec4 cur = texelFetch(text, coo ,0);
    return  slide(cur,tar,sli);
}

vec3 opU(in vec3 a, in vec3 b)
{
    return a.x < b.x ? a : b; 
}

float tri(float x)
{
    return min(fract(x) * 2., 2. - 2. * fract(x));
}

vec3 hash31( float n )
{
    return fract(sin(vec3(n,n+1.0,n+2.0))*vec3(43758.5453123,22578.1459123,19642.3490423));
}

mat2 rotate(float ang)
{
    return mat2(cos(ang), sin(ang),-sin(ang), cos(ang));
}

float smax( in float a, in float b, in float s ){
    float h = clamp( 0.5 + 0.5*(a-b)/s, 0.0, 1.0 );
    return mix(b, a, h) + h*(1.0-h)*s;
}
vec2 to_polar(vec3 U)
{
   return fract(vec2(
        atan(U.z, U.x) / 2.,
        atan(U.y, length(U.xz))
    ) / PI + .5);
}

vec3 to_cartesian(vec2 uv)
{
    uv = vec2(2, 1) * (uv - .5) * PI;
    return vec3(cos(uv.x), 0, sin(uv.x)) * cos(uv.y) + vec3(0, sin(uv.y), 0);
}

float sdCapsule( vec3 p, vec3 a, vec3 b, float r)
{
  vec3 pa = p - a, ba = b - a;
  float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
  return length( pa - ba*h ) - r;
}

vec2 polar(vec3 U)
{
   return fract(vec2(
        atan(U.z, U.x) / 2.,
        atan(U.y, length(U.xz))
    ) / PI + .5);
}
//==========================================================================================
//MIDI
//==========================================================================================
#define HEIGHT_CH_BLOCK 5

float getCCval (int CC, int channel,sampler2D midi_data)
{
    //add 1 to the channel to make it start from 1
    channel += 1;
    int   yco = (channel)*HEIGHT_CH_BLOCK-2;
    ivec2 coo = ivec2(CC, yco);
    float ccv = texelFetch(midi_data, coo,0).x;
    return ccv;
}

//Lighting Utils
float fresnel(float bias, float scale, float power, vec3 I, vec3 N)
{
    return bias + scale * pow(1.0 + dot(I, N), power);
}

//https://www.shadertoy.com/view/7lsBR8
// Very lame sky with sun-like object. Basically, this was a quick hack to emulate
// the "Forest Blurred" cube map... It needs work. :)
vec3 getSky(vec3 rd, vec3 ld){

    float lgt = max(dot(rd, ld), 0.);
    vec3 sky = mix(vec3(.1, .05, .04), vec3(.95, .97, 1)*3., clamp((rd.y + .5),0., 1.))/1.5;
    sky = mix(sky, vec3(8), pow(lgt, 8.));
    return min(sky*vec3(1, 1.1, 1.3), 1.);
}
//==========================================================================================
//ONESHADE code
//==========================================================================================
//https://www.shadertoy.com/view/slXXD4
// Definite integral of the light volume over the entire view of the ray (0 to ∞)
// I computed the limit based on the fact that atan() approaches π/2 (90 degrees)
// as x (the ratio of y over x) goes to infinity
float integrateLightFullView(in vec3 ro, in vec3 rd, in float k, in float d) 
{
    float a = dot(rd, rd);
    float b = dot(ro, rd);
    float c = dot(ro, ro) + d * d;
    float h = sqrt(a * c - b * b);
    return d * d * k * (RHO - atan(b / h)) / h;
}
//==========================================================================================
//NR4 code
//==========================================================================================

const vec3 c = vec3(1.,0.,-1.);
const float pi = acos(-1.);
// Determine zeros of k.x*x^2+k.y*x+k.z
vec2 quadratic_zeros(vec3 k)
{
    if(k.x == 0.) return -k.z/k.y*c.xx; // is a rare edgecase
    float d = k.y*k.y-4.*k.x*k.z;
    if(d<0.) return vec2(1.e4);
    return (c.xz*sqrt(d)-k.y)/(2.*k.x);
}
// Determine zeros of k.x*x^3+k.y*x^2+k.z*x+k.w
vec3 cubic_zeros(vec4 k)
{
    if(k.x == 0.) return quadratic_zeros(k.yzw).xyy;
    
    // Depress
    vec3 ai = k.yzw/k.x;
    
    //discriminant and helpers
    float tau = ai.x/3., 
        p = ai.y-tau*ai.x, 
        q = -tau*(tau*tau+p)+ai.z,
        dis = q*q/4.+p*p*p/27.;
        
    //triple real root
    if(dis > 0.) {
        vec2 ki = -.5*q*c.xx+sqrt(dis)*c.xz, 
            ui = sign(ki)*pow(abs(ki), c.xx/3.);
        return vec3(ui.x+ui.y-tau);
    }
    
    //three distinct real roots
    float fac = sqrt(-4./3.*p), 
        arg = acos(-.5*q*sqrt(-27./p/p/p))/3.;
    return c.zxz*fac*cos(arg*c.xxx+c*pi/3.)-tau;
}
// Point on a spline
vec3 xspline3(vec3 x, float t, vec3 p0, vec3 p1, vec3 p2)
{
    return mix(mix(p0,p1,t),mix(p1,p2,t),t);
}
// Distance to a point on a spline
float dspline3(vec3 x, float t, vec3 p0, vec3 p1, vec3 p2)
{
    return length(x - xspline3(x, t, p0, p1, p2));
}
// spline parameter of the point with minimum distance on the spline and sdf
// Returns vec2(dmin, tmin).
vec2 dtspline3(vec3 x, vec3 p0, vec3 p1, vec3 p2)
{
    vec3 E = x-p0, F = p2-2.*p1+p0, G = p1-p0;
    E = clamp(cubic_zeros(vec4(dot(F,F), 3.*dot(G,F), 2.*dot(G,G)-dot(E,F), -dot(E,G))),0.,1.);
    F = vec3(dspline3(x,E.x,p0,p1,p2),dspline3(x,E.y,p0,p1,p2),dspline3(x,E.z,p0,p1,p2));
    return F.x < F.y && F.x < F.z
        ? vec2(F.x, E.x)
        : F.y < F.x && F.y < F.z
            ? vec2(F.y, E.y)
            : vec2(F.z, E.z);
}
// Normal in a point on a spline
vec3 nspline3(vec3 x, float t, vec3 p0, vec3 p1, vec3 p2)
{
    return normalize(mix(p1-p0, p2-p1, t));
}
//distance to spline with parameter t
float dist(vec3 p0,vec3 p1,vec3 p2,vec3 x,float t)
{
    t = clamp(t, 0., 1.);
    return length(x-pow(1.-t,2.)*p0-2.*(1.-t)*t*p1-t*t*p2);
}
float dbox3(vec3 x, vec3 b)
{
  b = abs(x) - b;
  return length(max(b,0.))
         + min(max(b.x,max(b.y,b.z)),0.);
}
// Compute an orthonormal system from a single vector in R^3
mat3 ortho(vec3 d)
{
    vec3 a = normalize(
        d.x != 0. 
            ? vec3(-d.y/d.x,1.,0.)
            : d.y != 0.
                ? vec3(1.,-d.x/d.y,0.)
                : vec3(1.,0.,-d.x/d.z)
    );
    return mat3(d, a, cross(d,a));
}

vec2 asphere(vec3 x, vec3 dir, float R)
{
    float a = dot(dir,dir),
        b = 2.*dot(x,dir),
        cc = dot(x,x)-R*R,
        dis = b*b-4.*a*cc;
    if(dis<0.) return vec2(1000.);
    vec2 dd = (c.xz*sqrt(dis)-b)/2./a;
    return vec2(min(dd.x, dd.y), max(dd.x, dd.y));
}
/*
===================================================================================
*/

float sdCappedCylinder( vec3 p, float h, float r )
{
  vec2 d = abs(vec2(length(p.xz),p.y)) - vec2(h,r);
  return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

float smin( float a, float b, float k , inout float h ) 
{
    h = clamp( 0.5+0.5*(b-a)/k, 0., 1. );
    return mix( b, a, h ) - k*h*(1.0-h);
}

vec3 ropePerhi(vec3 p, int rope_id, sampler2D text, float env, inout vec3 hit_point)
{
    float dist = 1000.;
    vec2 uv;
    for(int i = 0; i < ROPE_POINTS-2; i+=2)
    {
        vec3 pp1  = texelFetch(text,ivec2(i+0,rope_id),0).xyz;
        vec3 pp2  = texelFetch(text,ivec2(i+1,rope_id),0).xyz;
        vec3 pp3  = texelFetch(text,ivec2(i+2,rope_id),0).xyz;
        if(length(pp1-pp2) < 0.001) return vec3(100);//avoid mapping static objects
        vec2 b = dtspline3(p,pp1,pp2,pp3);
        float t = b.y;
        vec3 norm  = nspline3(p,t,pp1,pp2,pp3);
        vec3 point = xspline3(p,t,pp1,pp2,pp3);
        mat3 m = ortho(norm),
            mt = transpose(m);
        vec3 c_point = mt*(p-point);
        
        const float span = 1./4.;
        float bez_ind = floor(float(i)/2.);
        vec2 lwise =vec2(t*span+bez_ind*span,t);
        float l = clamp(pow(lwise.x,0.8),0.,1.);
        float w =smoothstep(0.4,0.,abs(l-pow(env,0.2)))*smoothstep(0.2,0.6,l)*0.4+0.2;
        
        //float d = dbox3(c_point, vec3(.1,w*0.4+0.11, w*0.4+0.19));
         float d = sdCappedCylinder(c_point,w,w);
         d = max(d,length(c_point)-w*0.9);
        if(d < dist) 
        {
            hit_point = c_point;
            uv = vec2(c_point.z*w+w*0.5,lwise.x);
        }
        dist =  min(dist,d); 
    }
    return vec3(dist, uv);
}

vec3 ropePong(vec3 p, int rope_id, sampler2D text, float env, inout vec3 hit_point)
{
    float dist = 1000.;
    vec2 uv;
    for(int i = 0; i < ROPE_POINTS-2; i+=2)
    {
        vec3 pp1  = texelFetch(text,ivec2(i+0,rope_id),0).xyz;
        vec3 pp2  = texelFetch(text,ivec2(i+1,rope_id),0).xyz;
        vec3 pp3  = texelFetch(text,ivec2(i+2,rope_id),0).xyz;
        if(length(pp1-pp2) < 0.001) return vec3(100);//avoid mapping static objects
        vec2 b = dtspline3(p,pp1,pp2,pp3);
        float t = b.y;
        vec3 norm  = nspline3(p,t,pp1,pp2,pp3);
        vec3 point = xspline3(p,t,pp1,pp2,pp3);
        mat3 m = ortho(norm),
            mt = transpose(m);
        vec3 c_point = mt*(p-point);
        
        const float span = 1./4.;
        float bez_ind = floor(float(i)/2.);
        vec2 lwise =vec2(t*span+bez_ind*span,t);
        env = pow(env,0.5);

        float w =smoothstep(0.5,0.,abs(lwise.x-env))*0.4*smoothstep(0.9,0.1,lwise.x)+0.2;//+0.1*smoothstep(0.2,0.4,lwise.x);
        
        //float d = dbox3(c_point, vec3(.1, w, w));
         float d = sdCappedCylinder(c_point,w,w);
         d = max(d,length(c_point)-w*0.9);
        if(d < dist) 
        {
            hit_point = c_point;
            uv = vec2(c_point.z*w+w*0.5,lwise.x);
        }
        float h;
        dist =  min(dist,d);  
    }
    return vec3(dist, uv);
}

vec3 ropeDrums(vec3 p, int rope_id, sampler2D text, float env, inout vec3 hit_point)
{
    float dist = 1000.;
    vec2 uv;
    for(int i = 0; i < ROPE_POINTS-2; i+=2)
    {
        vec3 pp1  = texelFetch(text,ivec2(i+0,rope_id),0).xyz;
        vec3 pp2  = texelFetch(text,ivec2(i+1,rope_id),0).xyz;
        vec3 pp3  = texelFetch(text,ivec2(i+2,rope_id),0).xyz;
        if(length(pp1-pp2) < 0.001) return vec3(100);//avoid mapping static objects
        vec2 b = dtspline3(p,pp1,pp2,pp3);
        float t = b.y;
        vec3 norm  = nspline3(p,t,pp1,pp2,pp3);
        vec3 point = xspline3(p,t,pp1,pp2,pp3);
        mat3 m = ortho(norm),
            mt = transpose(m);
        vec3 c_point = mt*(p-point);
        
        const float span = 1./4.;
        float bez_ind = floor(float(i)/2.);
        vec2 lwise =vec2(t*span+bez_ind*span,t);
        float l = pow(lwise.x,1.7);
        float w = smoothstep(0.4,0.,abs(l-pow(env,0.5)))*0.2+smoothstep(-0.3,0.95, l)*1.4*0.1+0.01;
        //float d = dbox3(c_point, vec3(.1, w, w*0.5));
        float d = sdCappedCylinder(c_point,w,w*0.5);
         d = max(d,length(c_point)-w*0.4);
        if(d < dist) 
        {
            hit_point = c_point;
            uv = vec2(c_point.z*w+w*0.5,lwise.x);
        }
        dist =  min(dist,d); 
    }
    return vec3(dist*1.05, uv);
}


vec3 rope_flower1(vec3 p, int rope_id, sampler2D text, float env, inout vec3 hit_point)
{
    float dist = 1000.;
    vec2 uv;
    for(int i = 0; i < ROPE_POINTS-2; i+=2)
    {
        vec3 pp1  = texelFetch(text,ivec2(i+0,rope_id),0).xyz;
        vec3 pp2  = texelFetch(text,ivec2(i+1,rope_id),0).xyz;
        vec3 pp3  = texelFetch(text,ivec2(i+2,rope_id),0).xyz;
        if(length(pp1-pp2) < 0.001) return vec3(100);//avoid mapping static objects
        vec2 b = dtspline3(p,pp1,pp2,pp3);
        float t = b.y;
        vec3 norm  = nspline3(p,t,pp1,pp2,pp3);
        vec3 point = xspline3(p,t,pp1,pp2,pp3);
        mat3 m = ortho(norm),
            mt = transpose(m);
        vec3 c_point = mt*(p-point);
        
        const float span = 1./4.;
        float bez_ind = floor(float(i)/2.);
        vec2 lwise =vec2(t*span+bez_ind*span,t);
        float l = pow(lwise.x,0.7);
        float w =smoothstep(0.4,0.,abs(l-(1.-pow(env,0.5))))*3.8*smoothstep(0.8,0.6, l)*0.12+0.1;
        //float l = smoothstep(0.2,0.7,lwise.x)+0.01;
        
        //float d = dbox3(c_point, vec3(.1, w, w));
        float d = sdCappedCylinder(c_point,w,w);
        d = max(d,length(c_point)-w*0.9);
        //float d = sdCapsule(p,c_point,c_point*1.1,w);
        if(d < dist) 
        {
            hit_point = c_point;
            uv = vec2(c_point.z*0.8+0.2,lwise.x);
        }
        float hh;
        dist =  min(dist,d);//smin(dist,d, 0.18);   
    }
    //uv not returned(add code)
    return vec3(dist, uv);
}



//==================================================================================================
//Animations
//==================================================================================================
vec3 getRO(ivec2 tex_coo, int song_sect, int chan, ivec4 cc, sampler2D feed, sampler2D midi)
{
    if(iFrame < 10) return vec3(0,0,10);
    vec3  cur = texelFetch(feed,tex_coo,0).xyz;
    float x = texelFetch(midi, ivec2(cc.x,chan*HEIGHT_CH_BLOCK+3),0).x,
          y = texelFetch(midi, ivec2(cc.y,chan*HEIGHT_CH_BLOCK+3),0).x,
          z = texelFetch(midi, ivec2(cc.z,chan*HEIGHT_CH_BLOCK+3),0).x,
       dist = texelFetch(midi, ivec2(cc.w,chan*HEIGHT_CH_BLOCK+3),0).x; 
    vec3 tar = vec3(0,0,10);
    float sli = 0.3;
    float revolutions = 4.;
    float ang = y * TAU * revolutions;
    float[8] drums;
    for(int i = 0; i < 8; i++)
    {
        drums[i] =pow(texelFetch(feed,ivec2(i,DRUMS_ENV_ROW),0).x*.1,.5)*0.1;
    }
    sli = 0.051*STREAM_SLIDE;
    tar = vec3(cos(-ang+drums[7]+drums[5]+drums[3]),z+drums[6]+drums[2],sin(ang-drums[0]-drums[4]))*dist*RO_DIST_MULT;
    return slide(cur,tar, sli);
}

bool is_flower_on(int song_sect) {return song_sect < 2 || song_sect == 8 || song_sect > 9;}
bool is_pong_on(int song_sect)   {return song_sect > 0 && song_sect < 4;}
bool is_perhi_on(int song_sect)  {return song_sect > 1 && all(notEqual(ivec3(song_sect),ivec3(3,11,12)));}
bool is_drums_on(int song_sect)  {return song_sect > 3 && all(notEqual(ivec2(song_sect),ivec2(6,11   )));}
bool is_bass_on(int song_sect)   {return any(equal(ivec3(song_sect),ivec3(7,8,9)));}

float vel_note_on(sampler2D midi, int channel, int pitch, inout bool on) 
{
    ivec2 p = ivec2(pitch, HEIGHT_CH_BLOCK * channel );
    float vel = texelFetch(midi, p, 0).x;
    float secs_since_note_on  = texelFetch(midi, p + ivec2(0, 1), 0).x;
    float secs_since_note_off = texelFetch(midi, p + ivec2(0, 2), 0).x;
    float env = 0.;
    on = secs_since_note_on < secs_since_note_off || secs_since_note_off < 0.;
    return on ? vel : 0.;
}

float vel_note_on_drums(sampler2D midi, int channel, int pitch1, int pitch2, inout bool on) 
{
    //I made a mistake in mapping drums velocity in the synth so the actual velocity
    // of the drum hit comes from a different midi note so I need one pitch to check 
    //note on/off and one pitch for the velocity
    ivec2 p1 = ivec2(pitch1, HEIGHT_CH_BLOCK * channel );
    ivec2 p2 = ivec2(pitch2, HEIGHT_CH_BLOCK * channel );
    float vel = texelFetch(midi, p2, 0).x;
    float secs_since_note_on  = texelFetch(midi, p2 + ivec2(0, 1), 0).x;
    float secs_since_note_off = texelFetch(midi, p2 + ivec2(0, 2), 0).x;
    float env = 0.;
    on = secs_since_note_on < secs_since_note_off || secs_since_note_off < 0.;
    return on ? vel : 0.;
}

int get_time_sector(float time, float width)
{
    int time_section = int(mod(time,width*TIME_UNIT)/TIME_UNIT);
    return time_section;
}

vec3 getEnvelope(int id, int chan, int row, float dur, float sli, sampler2D midi, sampler2D feedb)
{
        int pitch = int(getCCval(19,chan,midi)*127.+0.01) ;
        vec4  prev_data = iFrame < 10 ? vec4(0) : texelFetch(feedb,ivec2(id,row),0);
        float prev_env  = prev_data.x, 
                was_on  = prev_data.y;
        bool on = false;
        float vel = vel_note_on(midi, chan, pitch, on);
        float is_on = on ? 1. : 0.;
        bool trigger = on && was_on < 0.5;
        float env = 0.;
        if(trigger) 
        {
            env = vel;
        }
        else 
        {
            env = slide(prev_env, 0., 1./dur); 
        }
        return vec3(slide(prev_env,env,sli), is_on, pitch);
}

vec4 getEnvelopeSector(int id, int chan, int row, float dur, float sli, float width, sampler2D midi, sampler2D feedb)
{
        int pitch = int(getCCval(19,chan,midi)*127.+0.01) ;
        vec4  prev_data = iFrame < 10 ? vec4(0) : texelFetch(feedb,ivec2(id,row),0);
        float prev_env  = prev_data.x, 
                was_on  = prev_data.y,
                sector = prev_data.w;
        bool on = false;
        float vel = vel_note_on(midi, chan, pitch, on);
        float is_on = on ? 1. : 0.;
        bool trigger = on && was_on < 0.5;
        float env = 0.;
        if(trigger) 
        {
            env = vel;
            sector = float(get_time_sector(iTime, width));
        }
        else 
        {
            env = slide(prev_env, 0., 1./dur); 
        }
        return vec4(slide(prev_env,env,sli), is_on, pitch, sector);
}

vec3 getEnvelopeDrums(ivec2 tex_coo, int pitch1, int pitch2, int chan, float dur, float sli, sampler2D midi, sampler2D feedb)
{
        vec4  prev_data = iFrame < 10 ? vec4(0) : texelFetch(feedb,tex_coo,0);
        float prev_env = prev_data.x, 
              was_on   = prev_data.y;
        bool on = false;
        float vel = vel_note_on_drums(midi, chan, pitch1, pitch2, on);
        float is_on = on ? 1. : 0.;
        bool trigger = on && was_on < 0.5;
        float env =0.;
        if(trigger) 
        {
            env = vel;
        }
        else 
        {
            env = slide(prev_env, 0., 1./dur); 
        }
        return vec3(slide(prev_env,env,sli), is_on, pitch1);
}

vec3 getEnvelopeBass(ivec2 tex_coo, int pitch, int chan, float dur, float sli, sampler2D midi, sampler2D feedb)
{
        vec4  prev_data = iFrame < 10 ? vec4(0) : texelFetch(feedb,tex_coo,0);
        float prev_env  = prev_data.x, 
                was_on  = prev_data.y;
        bool on = false;
        float vel = vel_note_on(midi, chan, pitch, on);
        float is_on = on ? 1. : 0.;
        bool trigger = on && was_on < 0.5;
        float env = 0.;
        if(trigger) 
        {
            env = vel;
        }
        else 
        {
            env = slide(prev_env, 0., 1./dur); 
        }
        return vec3(slide(prev_env,env,sli), is_on, pitch);
}

int getSongSection(sampler2D midi)
{
    int section = int(getCCval(50,0,midi)*127.+0.1);
    return section;
}

vec3 getTrajectory(int id, sampler2D midi)
{
    const int CC_offset = 26;
    //xyz pos are on 3 CC for each object starting from cc 26
    vec3 pos = vec3(getCCval(CC_offset+0, id, midi), getCCval(CC_offset+1, id, midi), getCCval(CC_offset+2, id, midi));
    return pos;
}

vec3 getTrajectoryDrums(int id, int chan, sampler2D midi)
{
    const int CC_offset = 50;
    int CC_chan  = chan;
    //xyz pos are on 3 CC for each object
    //always on channel 16, starting from CC 50
    int CC_ind = CC_offset + id *3;
    //different CC ind: made a mistake in mapping velocity in synth
    //see vel_note_on_drums for explanation
    int CC_ind_vel = CC_offset + (NUM_DRUMS - id ) * 3;
    vec3 pos = vec3(getCCval(CC_ind+0, CC_chan, midi), 0.,
                    getCCval(CC_ind_vel+1, CC_chan, midi));
    return pos;
}

vec3 getTrajectoryBass(int id, int chan, sampler2D midi)
{
    const int CC_offset = 50;
    int CC_chan  = chan;
    //xyz pos are on 3 CC for each object
    //always on channel 15, starting from CC 50
    int CC_ind = CC_offset + id *3;
    //different CC ind: made a mistake in mapping velocity in synth
    //see vel_note_on_drums for explanation
    int CC_ind_vel = CC_offset +  id  * 3;
    vec3 pos = vec3(getCCval(CC_ind+0, CC_chan, midi),
                    getCCval(CC_ind+1, CC_chan, midi) ,
                    getCCval(CC_ind+2, CC_chan, midi)
                    );
    return pos;
}

vec3 getPosFlower(int id, float env, sampler2D midi)
{
    vec3 pos = vec3(0);
    //longitudinal offset for each object
    //data.x :  rhythmic, data.y : velocity, data.z : pitch
    vec3 data = getTrajectory(id, midi);
    //data.x = 1.-data.x;
    //data = pow(data,vec3(2.3));
    data.y = data.y*0.1+0.25;// pow(data.y,0.5);
    env = pow(env, 0.5);// smoothstep(0.051,0.3,env);
    //rhythmic position (within bar): how far from the center of space
    //closer to downbeat->closer to center
    //velocity position: longitude (rotating around center)
    //pitch position : latitude
    float longi = data.x*TAU+PI;
    float lati =  data.y*TAU*0.8;
    pos = vec3(cos(longi),cos(float(id)/5.*TAU)*0.7+pow(env,0.5)*1.4-0.2,sin(longi))*FLOWER_LENGTH;//*vec3(data.x,1,data.x)*4.+FULCRUM1;
    //pos = vec3(cos(data.z*TAU),float(id)/4.*2.-1.,sin(data.z*TAU))*3.;
    //pos = vec3( float(id),data.y*3.,0);

    return pos;
}

vec3 flower_sect_displ(int song_sect)
{
    return song_sect < 5 || song_sect > 10? vec3(0) : vec3(0,0,-LAST_SECT_DISPLACE);
}

vec3 getPosPong(int id, float env, sampler2D midi)
{
    vec3 pos = vec3(0);
    //longitudinal offset for each object
    //data.x :  rhythmic, data.y : velocity, data.z : pitch
    vec3 data = getTrajectory(id, midi);
    //data.x = 1.-data.x;
    //data = pow(data,vec3(0.5));
    //data.y = pow(data.y,0.5);
    env = smoothstep(0.051,0.6,env);
    float longi = data.x*TAU*0.1;
    float lati =  data.y*TAU*0.8;
    float offs = float(id-8);
    pos = vec3(data.x*0.2,0,offs+data.x*0.);
    pos.xz *= 0.25;
    //pos.z -= mod(float(id), 2.) > 0.5 ? 5. : 0.;
    vec3 spherical_pos = to_cartesian(pos.xz)*(STAR_RAD*1.-env*2.8+0.5); 
    //spherical_pos.y -= env*1.;
    spherical_pos.xz *= rotate(iTime);
    spherical_pos.zy *= rotate(iTime+0.2);

    return false ? pos : spherical_pos;
}
vec3 getPosPerhi(int id, vec4 data, sampler2D midi)
{
    float env = data.x, pitch = data.z/127.*2.-1., sect = data.w;
    const float perhi_time_width = 4.;
    float offs = float(id)/8.;
    float y = env > 0.01 ? pitch*1.1 : 0.;
    float z = sect;
    //vec3 s = vec3(cos(x*TAU+sect*0.01+iTime*0.5),y, sin(x*TAU+sect*0.01+iTime*0.5));
    vec3 s = vec3(sect/3.,float(id)/8.*2.-0.5+pitch*0.5,0);
    s = vec3(sect/3.*0.4+0.1,offs*2.8+0.08+pitch*0.1,0);
    //s.xy += vec2(0.2,0.1);
    //s *= vec3(0.05,1,1);
    //s.x += x*2.-1;
    return s;
}

vec3 perhi_sect_displ(int song_sect)
{
    int sect_check = song_sect < 4 ? 0 : song_sect > 6 ? 2 : 1; 
    return sect_check == 0 ? vec3(0,0,0) : 
           sect_check == 1 ? vec3(0) : 
           vec3(0,0,LAST_SECT_DISPLACE) ;
}

vec3 getPosDrums(int id, float env, sampler2D midi)
{
    float drums_time_width = DRUMS_TIME_WIDTH[id];
    float sect = env > 0.1 ? (get_time_sector(iTime,drums_time_width)+1.)/(drums_time_width+1.) : 0.;
    vec3 data = getTrajectoryDrums(id,DRUMS_CHAN,midi)*1.;
   //data.x = 1.-data.x;
    //data = pow(data,vec3(0.5));
    //data.y = pow(data.y,0.5);
    env = smoothstep(0.051,0.3,env);
    float offs = float(id)/8.;
    float ang = offs*TAU;//+sect*0.+data.x*0.1+iTime*0.12;
    float rad = STAR_RAD;
    vec3 pos = vec3(cos(offs),sin(offs),0)*STAR_RAD*(env*0.5+0.2);
    //pos.xz *=rotate(iTime*0.12);
    //pos.z -= mod(float(id), 2.) > 0.5 ? 5. : 0.;
    return pos;
}

vec3 getPosBass(int id, float env, sampler2D midi)
{
    vec3 data = getTrajectoryBass(id,BASS_CHAN,midi)*1.;
   //data.x = 1.-data.x;
    data = pow(data,vec3(0.5));
    //data.y = pow(data.y,0.5);
    env = smoothstep(0.051,0.3,env);
    float offs = float(id)/7.;
    float ang = offs*TAU+data.x*0.4+data.y;
    float rad = BASS_TUBE_RADIUS-pow(env,0.5)*1.8;
    vec3 pos = vec3(rad*cos(ang),rad*sin(ang),-2*data.z*2.);
    //pos.z -= mod(float(id), 2.) > 0.5 ? 5. : 0.;
    return pos;
}
vec3 animFlowerData(ivec2 tex_coo, ivec4 block, sampler2D midi, sampler2D text)
{
    int rope_points = block.x, sub_block = block.z, tot_block = block.y;
    float frope_points = float(rope_points);
    int block_id = int(float(tex_coo.y)/float(sub_block));
    const int chan_offset = 0;
    //return vec3(block_id == 4 ? 1 : 0);
    if(block_id == 0)
    {
        //free flowing positions, not grounded to the center
        //get pos on first row pixel and transmit it to following pixels 
        if (tex_coo.x == 0)
        {
            int id = tex_coo.y;
            float env = texelFetch(text, ivec2(id, FLOWER_ENV_ROW),0).x;
            return getPosFlower(id,env, midi);
        } 
        else  return iFrame < 10 ? vec3(1) : pix_stream(tex_coo,text,0.9*STREAM_SLIDE);
    }
    else if(block_id == 1)
    {
        //attach each stem to the center
        int id = tex_coo.y - sub_block*1;
        float stretch_ind = float(tex_coo.x)/frope_points;
        vec3 pos = texelFetch(text, tex_coo - ivec2(0,sub_block),0).xyz;
        pos = mix(vec3(0), pos , stretch_ind);
        return pos;
    }
    else if(block_id == 2)
    {
        //petals
        //stream from last position of previous stream
        int id = tex_coo.y - sub_block*2;
        float env = pow(texelFetch(text,ivec2(id, FLOWER_ENV_ROW),0).x,0.5);
        int pos_on_r = int(float(rope_points)*(1.-env*1.));
        if (tex_coo.x == 0) return texelFetch(text,ivec2(pos_on_r,id),0).xyz;
        else  return iFrame < 10 ? vec3(1) : pix_stream(tex_coo,text,0.6);
    }
    else if(block_id == 3)
    {
        //stored positions with stretch, stretch from direction vector instead of anchoring it to the center like 
        //previous stream
        float stretch_ind = float(tex_coo.x)/frope_points;
        vec3 pos  = texelFetch(text,tex_coo - ivec2(0,sub_block),0).xyz;
        vec3 fulc = vec3(0);
        vec3 dir = max(normalize(pos - fulc),vec3(0.01));//texelFetch(text,tex_coo + ivec2(1,NUM_FLOWER_PETALS),0).xyz;
        pos += dir*stretch_ind*0.;

        return  pos;
    }
    //this row is just to get the envelope
    else if(block_id == 4) 
    {
        int id = tex_coo.x;
        int chan = id + chan_offset;
        return getEnvelope(id, chan, FLOWER_ENV_ROW, 6., 0.9, midi, text);
    }
}

vec3 animPongData(ivec2 tex_coo, ivec4 block, sampler2D midi, sampler2D text)
{
    int rope_points = block.x, sub_block = block.z, tot_block = block.y, block_offset = block.w;
    float frope_points = float(rope_points);
    int block_id = int(float(tex_coo.y-block_offset)/float(sub_block));
    // pong is on midi chans 9-16
    const int chan_offset = 8;
    //return vec3(block_id == 1 ? 1 : 0);
    if(block_id == 0)
    {
        //free flowing positions, not grounded to the center
        //get pos on first row pixel and transmit it to following pixels 
        if (tex_coo.x == 0)
        {
            int id = tex_coo.y - block_offset ;
            float env = texelFetch(text, ivec2(id, PONG_ENV_ROW),0).x;
            vec3 tar = getPosPong(id+chan_offset,env, midi);
            vec3 cur = texelFetch(text,tex_coo,0).xyz;
            return slide(tar, cur, 0.5*STREAM_SLIDE);
        } 
        else  return iFrame < 10 ? vec3(1) : pix_stream(tex_coo,text,0.5*STREAM_SLIDE);
    }
    else if(block_id == 1)
    {
        //return vec3(0,0,0);
        int id = tex_coo.y - block_offset - sub_block;
        float stretch_ind = float(tex_coo.x)/frope_points;
        vec3 pos = texelFetch(text,tex_coo-ivec2(0,block.z),0).xyz;
        float x = float(id)/2.*2.-1.;
        pos = pos*(1.+stretch_ind*0.6);
        pos.xz *= rotate(PI);
        //pos.yz *= rotate(iTime*0.12+RHO);
        //pos.y += stretch_ind*5.;
        return pos;
    }
    else if(block_id == 2)
    {
        int id = tex_coo.x;
        int chan = id + chan_offset;
        return getEnvelope(id, chan, PONG_ENV_ROW, 3., 0.5, midi, text);
    } 

}

vec4 animPerhiData(ivec2 tex_coo, ivec4 block, sampler2D midi, sampler2D text)
{
    int rope_points = block.x, sub_block = block.z, tot_block = block.y, block_offset = block.w;
    float frope_points = float(rope_points);
    int block_id = int(float(tex_coo.y-block_offset)/float(sub_block));
    // perhi is on midi chans 5-12
    const int chan_offset = 0;
    //return vec3(block_id == 1 ? 1 : 0);
    if(block_id == 0)
    {
        //free flowing positions, not grounded to the center
        //get pos on first row pixel and transmit it to following pixels 
        if (tex_coo.x == 0)
        {
            int id = tex_coo.y - block_offset;
            //get env pitch and time sector
            vec4 data = texelFetch(text, ivec2(id, PERHI_ENV_ROW),0);
            vec3 cur = texelFetch(text,tex_coo,0).xyz;
            vec3 tar = getPosPerhi(id+chan_offset,data, midi);
            //add a bit of rotation so that the bezier is always 
            //moving and is drawn correclty (if still it does not trace correctly)
            tar.xz += tar.xz*rotate(iTime*2.083)*0.01;
            return vec4(slide(tar, cur, 0.5),0);
        } 
        else  return iFrame < 10 ? vec4(1) :  pix_stream4(tex_coo,text,0.8*STREAM_SLIDE);
    }
    else if(block_id == 1)
    {
        //return vec3(0,0,0);
        int id = tex_coo.y - block_offset - sub_block;
        float stretch_ind = float(tex_coo.x)/(frope_points-1.);
        vec3 pos = texelFetch(text,tex_coo-ivec2(0,block.z),0).xyz;
        pos = to_cartesian(pos.xy)*STAR_RAD*0.8;
        float x = float(id)/2.*2.-1.;
        pos = mix(pos,pos*vec3(0.),stretch_ind);
        return vec4(pos,0);
    }
    else if(block_id == 2)
    {
        int id = tex_coo.x;
        int chan = id + chan_offset ;
        float dur = 1.8, width =4.;
        return getEnvelopeSector(id, chan, PERHI_ENV_ROW, dur, 0.5, width, midi, text);
    } 
}

vec3 animDrumsData(ivec2 tex_coo, ivec4 block, sampler2D midi, sampler2D text)
{
    int rope_points = block.x, sub_block = block.z, tot_block = block.y, block_offset = block.w;
    float frope_points = float(rope_points);
    int block_id = int(float(tex_coo.y-block_offset)/float(sub_block));
    // drums is on midi chan 16
    //return vec3(block_id == 1 ? 1 : 0);
    if(block_id == 0)
    {
        //free flowing positions, not grounded to the center
        //get pos on first row pixel and transmit it to following pixels 
        if (tex_coo.x == 0)
        {
            int id = tex_coo.y - block_offset ;
            //get env pitch and time sector
            float env = texelFetch(text, ivec2(id, DRUMS_ENV_ROW),0).x;
            vec3 cur = texelFetch(text,tex_coo,0).xyz;
            vec3 tar = getPosDrums(id,env, midi);
            return slide(tar, cur, 0.5*STREAM_SLIDE);
        } 
        else  return iFrame < 10 ? vec3(1) :  pix_stream(tex_coo,text,0.5*STREAM_SLIDE);
    }
    else if(block_id == 1)
    {
        //return vec3(0,0,0);
        int id = tex_coo.y - block_offset - sub_block;
        float stretch_ind = float(tex_coo.x)/frope_points;
        vec3 pos = texelFetch(text,tex_coo-ivec2(0,block.z),0).xyz;
        float stretch = 1.;
        float x = float(id)/2.*2.-1.;
        float env = tri(texelFetch(text,ivec2(id,DRUMS_ENV_ROW),0).x );
        //pos = mix(pos,vec3(x,0,0),stretch_ind);
        pos *= 1.+stretch_ind*DRUMS_LENGTH;
        pos.z += stretch_ind*1.2;
        int song_sect  = getSongSection(midi);
        pos.z +=  song_sect > 5 ? flower_sect_displ(song_sect).z : 0.; 
        pos.z += env*10.-10.;
        //pos.xz *= rotate(RHO);       
        return pos;
    }
    else if(block_id == 2)
    {
        
        int pitch1 = tex_coo.x + DRUMS_PITCH_OFFSET;
        int pitch2 = (NUM_DRUMS- tex_coo.x) + DRUMS_PITCH_OFFSET;
        int chan = DRUMS_CHAN ;
        return getEnvelopeDrums(tex_coo, pitch1, pitch2-1, chan, 5., 0.6, midi, text);
    } 
}

vec3 animBassData(ivec2 tex_coo, ivec4 block, sampler2D midi, sampler2D text)
{
    int rope_points = block.x, sub_block = block.z, tot_block = block.y, block_offset = block.w;
    float frope_points = float(rope_points);
    int block_id = int(float(tex_coo.y-block_offset)/float(sub_block));
    // bass is on midi chan 15
    //return vec3(block_id == 1 ? 1 : 0);
    if(block_id == 0)
    {
        //free flowing positions, not grounded to the center
        //get pos on first row pixel and transmit it to following pixels 
        if (tex_coo.x == 0)
        {
            int id = tex_coo.y - block_offset ;
            //get env pitch and time sector
            float env = texelFetch(text, ivec2(id, BASS_ENV_ROW),0).x;
            vec3 cur = texelFetch(text,tex_coo,0).xyz;
            vec3 tar = getPosBass(id,env, midi);
            return slide(tar, cur, 0.5*STREAM_SLIDE);
        } 
        else  return iFrame < 10 ? vec3(1) :  pix_stream(tex_coo,text,0.5*STREAM_SLIDE);
    }
    else if(block_id == 1)
    {
        //return vec3(0,0,0);
        int id = tex_coo.y - block_offset - sub_block;
        float stretch_ind = float(tex_coo.x)/frope_points;
        vec3 pos = texelFetch(text,tex_coo-ivec2(0,block.z),0).xyz;
        float stretch = 1.;
        float x = float(id)/2.*2.-1.;
        //pos = mix(pos,vec3(x,0,0),stretch_ind);
        pos *= 1.+stretch_ind*BASS_LENGTH;
        pos.z += stretch_ind*1.;
        return pos;
    }
    else if(block_id == 2)
    {
        
        int pitch = tex_coo.x + BASS_PITCH_OFFSET;
        int chan = BASS_CHAN ;
        return getEnvelopeBass(tex_coo, pitch, chan, 5., 0.5, midi, text);
    } 
}
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

vec3 plane(vec3 ro, vec3 rd, vec3 p, vec3 norm)
{
    float t = max(0., dot(p-ro,norm)/dot(rd,norm));
    return ro+rd*t;
}

vec3 blinn_phong(vec3 p, vec3 rd, vec3 light, vec3 norm,  vec3 col_diffuse)
{
    vec3 col_ambient = vec3(0.8588, 0.8196, 0.098);
    vec3 col_specular = vec3(0.3686, 0.3725, 0.3059);
    return  col_ambient + 
            col_diffuse * max(dot(normalize(light-p),norm),0.)+ 
            col_specular * pow(max(dot(reflect(normalize(light-p),norm),rd),0.),2.);

}

vec3 get_light2(HitInfo hit, vec3 rd, vec3 l_pos, vec3 l_col, sampler2D TXT)
{
    vec3 norm= hit.nor, p= hit.pos;
    vec3 diff = blinn_phong(p,rd,SunPos,norm,SunCol)*0.1;
    vec3 r = reflect(rd,norm);
    float l = length(p);
    float fr = clamp(1. - dot(norm,-rd), 0.,1.);
    float bodyR = 1.5;
    float cone = cos(atan(bodyR / l));

    float specshad = 1. - smoothstep(-0.1,0.1, dot(r,normalize(-p)) -cone);
    float specshad2 = 1. - smoothstep(-0.3,0.3, dot(r,normalize(-p)) -cone);

    //backlight / fake SSS
    vec3 col = diff;//+vec3(1,1,0.4) * pow(env,1.)* 1.1/(10.-exp(hit.x*.01)); 
    //col += l_col * mix(vec3(0.2,0.5,1.) /  2., vec3(1.,0.9,0.8).bgr,specshad2 * pow(fr,0.8))*0.8;

    //fake AO from center 
    col *= vec3(pow(mix(0.5+0.5*dot(norm,normalize(-p)), 1., smoothstep(0.8,0.9,l)), 0.5));

    //slight AO / diffuse bleeding
    col *= mix(SunCol, vec3(1), smoothstep(0.5,0.9,l));

    vec3 c = col;
    //yellow tips
    //col = mix(c.bbb * vec3(0.5,1.,0.5), col *pow(hit.env,1.)*0.1, smoothstep(0.65,1.,l));
        
    //envelope illumination
    float env_ill = smoothstep(0.1,0.,abs(l-pow(hit.env,1.5)));
    col += vec3(1,1,0.4) *(0.1/(1.-exp(-hit.dist))* pow(hit.env,1.)*3.1*env_ill);
    //col += vec3(1,1,0.4) * pow(hit.env,2.)*0.1;
    //specular highlight
    col +=  smoothstep(0.4,0.7, dot(r, normalize(vec3(1)))) * fr;

    if(hit.id.x == 3) hit.env = 1.- hit.env;
    //float stripe = pow(smoothstep(0.1,0.0,abs(pow(hit.env,0.5)-(1.-hit.uv.y))),0.5);
    //col *=vec3(stripe *0.5+0.8);//diff*pow(hit.env,1.5)*1.+0.1;//max(col, vec3(0));
    vec3 txt = texture(TXT,hit.uv_transorm).xyz;
    float bump = clamp(dot(txt,txt),0.,1.);
    col = mix(col,txt,0.6);//*2.2 *light_intensity*0.5;
    //mist
    //col += vec3(mix(vec3(0.0275, 0.4784, 0.9922), vec3(0), exp(-hit.dist /1005.)));

    return col  ;
} 
