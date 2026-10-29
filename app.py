"""Entry point for Python, PyInstaller and the Windows embedded runtime."""
import os
from pathlib import Path
import sys
import tempfile

os.environ.setdefault('MPLCONFIGDIR',str(Path(tempfile.gettempdir())/'elspin-matplotlib'))

if __name__=='__main__':
    try:
        if '--self-test' in sys.argv:
            from electrospin.selfcheck import main as selfcheck
            selfcheck()
            raise SystemExit(0)
        from electrospin.desktop import main
        raise SystemExit(main())
    except Exception:
        import traceback
        log=Path(tempfile.gettempdir())/'ElSpin-startup-error.log'
        log.write_text(traceback.format_exc(),encoding='utf-8')
        if sys.platform=='win32':
            import ctypes
            ctypes.windll.user32.MessageBoxW(None,'Не удалось запустить ElSpin. Подробности: '+str(log),'ElSpin',0x10)
        raise
