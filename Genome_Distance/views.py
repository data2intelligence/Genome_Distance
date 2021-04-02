import os

from django.http.response import FileResponse
from django.shortcuts import render
from Genome_Distance.settings import data_path

def send_file(request, f):
    f = os.path.join(data_path, f)
    
    if os.path.exists(f):
        return FileResponse(open(f, 'rb'))
    else:
        return render(request, 'error.html', {'message': 'File does not exist %s' % os.path.basename(f)})

def index(request):
    return render(request, 'index.html')

def contact(request):
    return render(request, 'contact.html')

def download(request):
    return render(request, 'download.html')
