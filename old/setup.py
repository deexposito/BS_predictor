from setuptools import setup, find_packages

setup(
    name='BS_predictor',
    version='1.1.0',
    description = "Machine Learning-Based Binding Site Prediction for Protein Structures",
    
    author="Denis Expósito Navarro, Marta Alonso Caubilla, Yolanda Andrés López",
    author_email="deexna4@gmail.com, marta.alonso07@estudiant.upf.edu, yolanda.andres01@estudiant.upf.edu",
    
    packages=find_packages(),
    
    package_data = {
        "BS_predictor": ["data/*.csv", "data/*.py", "data/train_pdbs/*", "data/*.txt"],
    },
    include_package_data=True,
    
    entry_points={
        'console_scripts': [
            'BS_predictor=BS_predictor:main'
        ]
    },    
    scripts=['BS_predictor.py'])

