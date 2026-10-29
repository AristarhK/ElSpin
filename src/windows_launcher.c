/* Native Unicode GUI launcher for the official CPython embedded runtime.
 * No shell, system Python, PATH mutation, installer, or network is needed.
 */
#ifndef UNICODE
#define UNICODE
#endif
#ifndef _UNICODE
#define _UNICODE
#endif
#include <windows.h>
#include <wchar.h>

int WINAPI wWinMain(HINSTANCE instance,HINSTANCE previous,PWSTR command,int show) {
    (void)instance; (void)previous; (void)command; (void)show;
    static wchar_t base[32768], python[32768], args[65536];
    DWORD size=GetModuleFileNameW(NULL,base,32768);
    if (!size || size>=32700) return 1;
    wchar_t *last=wcsrchr(base,L'\\');
    if (!last) return 1;
    *last=L'\0';
    if (swprintf(python,32768,L"%ls\\_runtime\\pythonw.exe",base)<0) return 1;
    if (swprintf(args,65536,L"\"%ls\" -X utf8 \"%ls\\app.py\"",python,base)<0) return 1;
    STARTUPINFOW startup={0}; startup.cb=sizeof(startup);
    PROCESS_INFORMATION process={0};
    if (!CreateProcessW(python,args,NULL,NULL,FALSE,0,NULL,base,&startup,&process)) {
        MessageBoxW(NULL,L"Не удалось запустить встроенный Python. Распакуйте весь архив, включая папку _runtime.",L"ElSpin",MB_OK|MB_ICONERROR);
        return 1;
    }
    CloseHandle(process.hThread);WaitForSingleObject(process.hProcess,INFINITE);
    DWORD exit_code=1;GetExitCodeProcess(process.hProcess,&exit_code);CloseHandle(process.hProcess);
    return (int)exit_code;
}
