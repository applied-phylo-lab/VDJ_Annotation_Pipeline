import asyncio
import sys
import shutil
import os
import csv
import time


async def run_script(script_path: str, args: list[str], job_id: int) -> dict:
    time.sleep(2)
    proc = await asyncio.create_subprocess_exec(
        sys.executable, script_path, *args,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
    )
    
    stdout, stderr = await proc.communicate()
    
    time.sleep(2)
    return {
        "args": args,
        "returncode": proc.returncode,
        "stdout": stdout.decode().strip(),
        "stderr": stderr.decode().strip(),
    }

async def run_pool(script_path: str, all_args: list[list[str]], max_concurrent):
    queue = list(enumerate(all_args))         # [(id, args), ...]
    pending: set[asyncio.Task] = set()

    def launch_next():
        if queue:
            job_id, args = queue.pop(0)
            print(f"  → Starting job [{job_id}] with args: {args}")
            task = asyncio.create_task(
                run_script(script_path, args, job_id),
                name=f"job-{job_id}"
            )
            pending.add(task)

    # Fill the pool initially
    for _ in range(min(max_concurrent, len(all_args))):
        time.sleep(2)
        launch_next()

    while pending:
        done, _ = await asyncio.wait(pending, return_when=asyncio.FIRST_COMPLETED)

        for task in done:
            pending.discard(task)
            result = task.result()

            #if result["stdout"]:
            #    print(f"    stdout: {result['stdout'][-300:]}")
            #if result["args"]:
            #    print(result["args"])
            
            time.sleep(2)
            launch_next()  # immediately fill the vacant slot

    print("\nAll jobs finished.")

hep_thresh_range=sys.argv[1] # ex) -h:1-3
non_thresh_range=sys.argv[2] # ex) -n:2-6
thresh_exclude=sys.argv[3] # ex) -e:3-6,3-5 (-all to use all threshes)
thresh_type=sys.argv[4] # ex) -base:IGH
both_same=sys.argv[5] # ex) -both_same:1-9 (-none to not use)
motif=sys.argv[6] # ex) motif:1-2:prop (-none to not use)
size=sys.argv[7] # ex) -size:2 (-none to not use)
output=sys.argv[8] # ex) -outdel:all_data/all_birds_igh (-none to not save to a folder)
data_location=sys.argv[9] # ex) -csv (provide gene_list.csv in input_data)  |  -fold : path to IgDetective data folder with (order --> species --> haplotype) subdir structure
filtering=sys.argv[10] # ex) -f:bStrDea1_pri (-n for none)

if output!="-none":
    if str(output).startswith("-outdel") and os.path.isdir(output.split(":")[1].split("/")[0]):
        shutil.rmtree(output.split(":")[1].split("/")[0])

exclude_threshes=[]
if thresh_exclude!="-all":
    exclude_threshes = thresh_exclude.split(":")[1].split(",")

async def main():
    script = "d_gene_search.py"

    jobs=[]
    hep_range=[int(hep_thresh_range.split(":")[1].split("-")[0]),int(hep_thresh_range.split(":")[1].split("-")[1])]
    for h in range((hep_range[-1]-hep_range[0])-1):
        hep_range.append((hep_range[0]+h+1))

    non_range=[int(non_thresh_range.split(":")[1].split("-")[0]),int(non_thresh_range.split(":")[1].split("-")[1])]
    for n in range((non_range[-1]-non_range[0])-1):
        non_range.append((non_range[0]+n+1))


    for hep in hep_range:
        for non in non_range:
            thresh_arg=thresh_type.replace(":",":"+str(hep)+"-"+str(non)+":")
            if thresh_arg.split(":")[1] not in exclude_threshes:
                jobs.append([thresh_arg,both_same,motif,size,output+"_"+str(hep)+"-"+str(non),data_location,filtering])

    print(f"Running {len(jobs)} instances of {script}, up to 16 at a time...\n")
    await run_pool(script, jobs, max_concurrent=16)

asyncio.run(main())