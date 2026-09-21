$ErrorActionPreference = 'Continue'
Write-Output '=== WINDOWS ==='
Get-CimInstance Win32_OperatingSystem | Select-Object Caption, Version, BuildNumber
Write-Output '=== RAM (GiB) ==='
[math]::Round((Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory / 1GB, 1)
Write-Output '=== FREE DISK ==='
Get-PSDrive -PSProvider FileSystem | Select-Object Name, @{Name='FreeGiB';Expression={[math]::Round($_.Free/1GB,1)}}
Write-Output '=== WSL ==='
wsl --status
wsl --list --verbose
Write-Output '=== NVIDIA ==='
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv
