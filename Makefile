.PHONY: setup pipeline dashboard

setup:
	pip install -r requirements.txt

pipeline:
	python load_data.py
	python analysis.py

dashboard:
	streamlit run dashboard.py --server.headless true --server.port 8501
