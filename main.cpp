#include <string>
#include <iostream>

#include <stdint.h>

#include "CLI/App.hpp"
#include "CLI/Formatter.hpp"
#include "CLI/Config.hpp"

#define WIN32_LEAN_AND_MEAN
#define WIN32_EXTRA_LEAN
#include <Windows.h>

#include "revision_submission.h"

#define WIDTH 1920
#define HEIGHT 1080

typedef GLuint (APIENTRYP PFNGLCREATESHADERPROGRAMVPROC) (GLenum type, GLsizei count, const GLchar *const*strings);

#define ATOM_STATIC 0xc019

typedef struct {
	BYTE dmDeviceName[CCHDEVICENAME];
	WORD dmSpecVersion;
	WORD dmDriverVersion;
	WORD dmSize;
	WORD dmDriverExtra;
	DWORD dmFields;
	union {
		struct {
			short dmOrientation;
			short dmPaperSize;
			short dmPaperLength;
			short dmPaperWidth;
			short dmScale;
			short dmCopies;
			short dmDefaultSource;
			short dmPrintQuality;
		};
	
		struct {
			POINTL dmPosition;
			DWORD dmDisplayOrientation;
			DWORD dmDisplayFixedOutput;
		};
	};

	short dmColor;
	short dmDuplex;
	short dmYResolution;
	short dmTTOption;
	short dmCollate;
	BYTE dmFormName[CCHFORMNAME];
	WORD dmLogPixels;
	DWORD dmBitsPerPel;
	DWORD dmPelsWidth;
	DWORD dmPelsHeight;
	union {
		DWORD dmDisplayFlags;
		DWORD dmNup;
	};

	DWORD dmDisplayFrequency;
} HGDEVMODE;

#pragma data_seg("devmode")
static HGDEVMODE devmode = {
	{0}, 0, 0, sizeof(HGDEVMODE), 0, DM_PELSWIDTH | DM_PELSHEIGHT,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, {0}, 0, 0, WIDTH, HEIGHT, 0, 0
};

#pragma data_seg("pfd")
static const PIXELFORMATDESCRIPTOR pfd = {
	sizeof(PIXELFORMATDESCRIPTOR), 1,
	PFD_DRAW_TO_WINDOW|PFD_SUPPORT_OPENGL|PFD_DOUBLEBUFFER,
	PFD_TYPE_RGBA, 32, 0, 0, 0, 0, 0, 0, 8, 0,
	0, 0, 0, 0, 0, 32, 0, 0, PFD_MAIN_PLANE, 0, 0, 0, 0
};

int main(int argc, char **args)
{
	// Command line argument parser.
	float t_from = 0.,
		t_to = 100.;
	bool real_time = false,
		loop = false;
	int fps = 30;
	std::string output_dir = "output";
	
	CLI::App app{"Hyperborean sun renderer. Render this shiney."};
	app.add_option("-f,--from", t_from, "Time to start rendering.");
	app.add_option("-t,--to", t_to, "Time to stop rendering.");
	app.add_flag("-r,--real-time", real_time, "Play the video in real-time.");
	app.add_flag("-l,--loop", loop, "Loop between from-time and to-time.");
	app.add_option("--fps", fps, "Set the target fps for prerendering.");
	app.add_option("-o,--output", output_dir, "Save the output to this directory.")
		->check(CLI::ExistingDirectory);

	CLI11_PARSE(app, argc, args);

	// Create the window.
	ChangeDisplaySettingsA((DEVMODEA*)&devmode, CDS_FULLSCREEN);
	HWND hWnd = CreateWindowExA(0, (char*)ATOM_STATIC, 0, WS_POPUP|WS_VISIBLE|WS_MAXIMIZE, 0, 0, 0, 0, 0, 0, 0, 0);
	ShowCursor(0);

	HDC hDC = GetDC(hWnd);
	SetPixelFormat(hDC, ChoosePixelFormat(hDC, &pfd), &pfd);
	wglMakeCurrent(hDC, wglCreateContext(hDC));

	// Set up the shader.
	SETUP_REVISION_SUBMISSION(WIDTH, HEIGHT);
	SETUP_MUSIC;
	SETUP_SCREENSHOTS(WIDTH, HEIGHT);

#ifdef HAS_MIDI
	SETUP_MIDI;
#endif // HAS_MIDI

	if(real_time) {
		PLAY_MUSIC;
	}

	// Render loop.
	float iTime = t_from, last_time = -1.;
	for(int iFrame = 0;; ++iFrame) {
		if(real_time) {
			UPDATE_MUSIC_PLAYBACK_TIME;
			iTime = MUSIC_PLAYBACK_TIME;
		} else {
			iTime += 1./(float)fps;
		}

#ifdef HAS_MIDI
		UPDATE_MIDI_TEXTURE(iTime, last_time);
#endif // HAS_MIDI

		DRAW_REVISION_SUBMISSION(WIDTH, HEIGHT, iTime, iFrame);

		SwapBuffers(hDC);

		if(!real_time) {
			SCREENSHOT(WIDTH, HEIGHT, iFrame);
		}

		last_time = iTime;

		// Message loop
		MSG msg = {0};
		while (PeekMessageA(&msg, NULL, 0, 0, PM_REMOVE))
		{
			if(msg.message == WM_KEYDOWN && msg.wParam == VK_ESCAPE)
			{
				return -1;
			}
			
			DispatchMessageA(&msg);
		}
	}

	return 0;
}