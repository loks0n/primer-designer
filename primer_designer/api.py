#!/usr/bin/env python
"""
FastAPI interface for primer-designer
"""

from typing import List, Dict, Any
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
import uvicorn

from primer_designer.core import design_primers
from primer_designer.sequence import validate_dna
from primer_designer.cli import parse_mutations

app = FastAPI(
    title="Primer Designer API",
    description="API for designing primers for site-directed mutagenesis",
    version="0.1.0",
)


class PrimerDesignRequest(BaseModel):
    sequence: str
    mutations: List[str]
    min_flank_length: int = 15
    max_flank_length: int = 18
    min_tm: float = 65.0
    max_tm: float = 75.0
    optimal_tm: float = 70.0
    max_tm_diff: float = 5.0
    min_gc: float = 35.0
    max_gc: float = 70.0
    optimal_gc: float = 50.0


@app.post("/primers", response_model=List[Dict[str, Any]])
async def design_primers_endpoint(request: PrimerDesignRequest) -> List[Dict[str, Any]]:
    """
    Design primers for site-directed mutagenesis

    Parameters:
    - sequence: DNA sequence to use for primer design
    - mutations: List of mutations in format F10A (original AA, position, new AA)
    - min_flank_length: Minimum nucleotides flanking each side of mutation (default: 15)
    - max_flank_length: Maximum nucleotides flanking each side of mutation (default: 18)
    - min_tm: Minimum acceptable melting temperature in °C (default: 65.0)
    - max_tm: Maximum acceptable melting temperature in °C (default: 75.0)
    - optimal_tm: Target optimal melting temperature in °C (default: 70.0)
    - max_tm_diff: Maximum allowed difference between forward/reverse Tm (default: 5.0)
    - min_gc: Minimum acceptable GC content in percent (default: 35.0)
    - max_gc: Maximum acceptable GC content in percent (default: 70.0)
    - optimal_gc: Target optimal GC content in percent (default: 50.0)
    """
    try:
        # Validate DNA sequence
        sequence = validate_dna(request.sequence)

        # Parse mutations
        try:
            mutations = parse_mutations(request.mutations)
        except Exception as e:
            raise HTTPException(
                status_code=400, detail=f"Failed to parse mutations: {str(e)}"
            )

        # Design primers
        results = design_primers(
            sequence,
            mutations,
            min_flank_length=request.min_flank_length,
            max_flank_length=request.max_flank_length,
            min_tm=request.min_tm,
            max_tm=request.max_tm,
            optimal_tm=request.optimal_tm,
            max_tm_diff=request.max_tm_diff,
            min_gc=request.min_gc,
            max_gc=request.max_gc,
            optimal_gc=request.optimal_gc,
        )

        return results

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"An error occurred: {str(e)}")


if __name__ == "__main__":
    uvicorn.run("primer_designer.api:app", host="0.0.0.0", port=8000, reload=True)
