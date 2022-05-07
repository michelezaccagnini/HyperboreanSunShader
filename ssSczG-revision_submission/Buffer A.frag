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
    ivec2 blocks_dim = ivec2(9,
                            FLOWER_BLOCK.y+PONG_BLOCK.y+PERHI_BLOCK.y+PERLO_BLOCK.y+DRUMS_BLOCK.y+RO_BLOCK.y+BASS_BLOCK.y);                            
    if(any(greaterThanEqual(iCoo,blocks_dim))) discard;
    //if(iCoo.y == DRUMS_BLOCK_OFFSET) fragColor = vec4(0.5);
    else{
    int song_sect = getSongSection(MIDI);
    bool flower_on = is_flower_on(song_sect) && iCoo.y < PONG_BLOCK_OFFSET ;
    bool pong_on   = is_pong_on(song_sect) && iCoo.y >= PONG_BLOCK_OFFSET && iCoo.y < PERHI_BLOCK_OFFSET;
    bool perhi_on  = is_perhi_on(song_sect) && iCoo.y >= PERHI_BLOCK_OFFSET && iCoo.y < PERLO_BLOCK_OFFSET;
    bool drums_on  = is_drums_on(song_sect) && iCoo.y >= DRUMS_BLOCK_OFFSET && iCoo.y < RO_BLOCK_OFFSET;
    bool ro_on     = iCoo == ivec2(0,RO_COO.y);
    bool bass_on  = is_bass_on(song_sect) && iCoo.y >= BASS_BLOCK_OFFSET;
    if (ro_on)          fragColor.xyz =  true ? getRO(iCoo, song_sect, RO_CHAN, RO_CC, FEEDBACK, MIDI) : vec3(0,1,0);
    else if(flower_on)  fragColor.xyz =  true? animFlowerData(iCoo,FLOWER_BLOCK,MIDI,FEEDBACK):vec3(0.11);     
    else if(pong_on)    fragColor.xyz =  true? animPongData(iCoo,PONG_BLOCK,MIDI,FEEDBACK)  : vec3(1,0,0);//
    else if(perhi_on)   fragColor     =  true? animPerhiData(iCoo,PERHI_BLOCK,MIDI,FEEDBACK) : vec4(0,1,0,0);//
    else if(drums_on)   fragColor =  true ? animDrumsData(iCoo,DRUMS_BLOCK,MIDI,FEEDBACK) : vec4(0.451, 0.902, 0.902, 0.0);//
    else if(bass_on)    fragColor.xyz =  true ? animBassData(iCoo,BASS_BLOCK,MIDI,FEEDBACK) : vec3(0,0.2,1);//    
    }
}
#endif