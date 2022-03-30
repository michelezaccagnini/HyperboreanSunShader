//MIDI data parsing

#define FEEDBACK iChannel0
#define MIDI iChannel1
#define BUF_C iChannel2


//debug: visualize midi texture
#if 0
void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    fragColor = texelFetch(MIDI,ivec2(fragCoord),0);
}
#else

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    ivec2 iCoo = ivec2(fragCoord);
    ivec2 blocks_dim = ivec2(max(max(max(FLOWER_BLOCK.x,PONG_BLOCK.x),PERHI_BLOCK.x),PERLO_BLOCK.x),
                            FLOWER_BLOCK.y+PONG_BLOCK.y+PERHI_BLOCK.y+PERLO_BLOCK.y+DRUMS_BLOCK.y+RO_BLOCK.y);                            
    if(any(greaterThanEqual(iCoo,blocks_dim))) discard;
    //if(iCoo.y == DRUMS_BLOCK_OFFSET) fragColor = vec4(0.5);
    else{
    int song_sect = getSongSection(MIDI);
    bool flower_on = (song_sect < 2 || song_sect > 4) && iCoo.y < FLOWER_BLOCK.y;
    bool pong_on   = (song_sect > 0 && song_sect < 4 && iCoo.y >= PONG_BLOCK_OFFSET && iCoo.y < PERHI_BLOCK_OFFSET);
    bool perhi_on  = (song_sect > 1 && song_sect < 5 && iCoo.y >= PERHI_BLOCK_OFFSET && iCoo.y < PERLO_BLOCK_OFFSET);
    bool perlo_on  = song_sect > 3 && song_sect < 5 && iCoo.y >= PERLO_BLOCK_OFFSET && iCoo.y < DRUMS_BLOCK_OFFSET;
    bool drums_on  = song_sect > 3 && iCoo.y >= DRUMS_BLOCK_OFFSET && iCoo.y < RO_BLOCK_OFFSET;
    bool ro_on     = iCoo == ivec2(0,blocks_dim.y-1);
    if (ro_on)   fragColor.xyz = getRO(iCoo, song_sect, RO_CHAN, RO_CC, FEEDBACK, MIDI);
    else if(flower_on) fragColor.xyz =       true ? animFlowerData(iCoo,FLOWER_BLOCK,MIDI,FEEDBACK):vec3(1);     
    else if(pong_on)   fragColor.xyz =  true ? animPongData(iCoo,PONG_BLOCK,MIDI,FEEDBACK)  : vec3(1,0,0);//
    else if(perhi_on)  fragColor =  true ? animPerhiData(iCoo,PERHI_BLOCK,MIDI,FEEDBACK) : vec4(0,1,0,0);//
    //else if(perlo_on)  fragColor.xyz =  true ? animPerloData(iCoo,PERLO_BLOCK,MIDI,FEEDBACK) : vec3(1,1,0);//
    else if(drums_on)  fragColor.xyz =  true ? animDrumsData(iCoo,DRUMS_BLOCK,MIDI,FEEDBACK) : vec3(0,1,1);//
    
    }
}
#endif