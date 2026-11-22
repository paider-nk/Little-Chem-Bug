from pathlib import Path
import json
import tempfile
import unittest
from unittest.mock import patch

from gaussian_core import InputError, Settings, generate
from gaussian_workflow import run_jobs


class MascotProgressTests(unittest.TestCase):
    def test_generator_stages_follow_actual_preprocessing(self):
        settings = Settings(conformers=1, preprocess_method="none")
        stages = []
        result = generate("CCO", "smiles", settings, on_stage=stages.append)
        self.assertEqual(stages, ["解析结构与计算参数", "生成三维构象", "检查起始几何并组装输入文件"])
        self.assertEqual(result.text, generate("CCO", "smiles", settings).text)

    def test_invalid_input_does_not_report_geometry(self):
        stages = []
        with self.assertRaises(InputError):
            generate("!", "smiles", on_stage=stages.append)
        self.assertEqual(stages, ["解析结构与计算参数"])

    def test_batch_counts_completed_jobs_and_abnormal_exits(self):
        updates = []
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            inputs = [root / "first.gjf", root / "second.gjf"]
            for path in inputs:
                path.touch()

            def run(_, path, *_args):
                return {"input": str(path), "status": "incomplete" if path == inputs[0] else "normal"}

            output = root / "out"
            with patch("gaussian_workflow.locate_gaussian", return_value="g16.exe"), \
                 patch("gaussian_workflow._run_one", side_effect=run):
                results = run_jobs(inputs, output, cores=1, slots=1, on_progress=lambda *counts: updates.append(counts))
            self.assertEqual(updates[0], (0, 2, 0))
            self.assertEqual(updates[1][:2], (1, 2))
            self.assertEqual(updates[2], (2, 2, 1))
            self.assertEqual(results, json.loads((output / "run_report.json").read_text(encoding="utf-8")))


if __name__ == "__main__":
    unittest.main()
