PC Controller DOS Compatibility
The PC controller software for the Oncore family of GPS receivers uses a 
set of utility programs to handle the communications port setup and 
takedown.  These programs store the communications port settings when the 
controller program is started and then restore the original communications 
port settings when the controller program is terminated.  These utility 
programs only work with DOS versions 3, 4, or 5.  In order for them to work
properly when using DOS 6.x, the ‘SETVER’ command should be used.  This is
explained in this application note.
1)	Add the following line to the ‘CONFIG.SYS’ file:
		DEVICEHIGH /L:1,12048 = C:\DOS\SETVER.EXE
2)	Execute the following commands once at the DOS prompt (do not add 
to the ‘AUTOEXEC.BAT’ or ‘CONFIG.SYS’ files):
		C:\> SETVER MARKNET.EXE 4.00
		C:\> SETVER RELNET.EXE 4.00
		C:\> SETVER GPSCOMMS.EXE 4.00
3)	Reboot the computer.
The commands in step #2 add entries into the ‘SETVER’ database so that when
‘SETVER’ is loaded (in the ‘CONFIG.SYS’ file), the correct versions of DOS 
will be used when those programs are run.
