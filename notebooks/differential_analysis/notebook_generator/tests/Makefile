# ---- Config ----
CONFIG=configs/test.yaml
NOTEBOOK=notebooks/TestNotebook.ipynb
BLOCK=blocks/test/01_hello.json
SCRIPT_DIR=run

# ---- Targets ----

all: test

edit-notebook:
	@echo "✏️  Modifying notebook for reverse sync test"
	pwd
	sed -i "s/Goal/Way/" $(NOTEBOOK)

reverse-sync:
	@echo "🔄 Syncing notebook → block"
	../$(SCRIPT_DIR)/sync.sh $(NOTEBOOK) --update

test:
	@echo "🧪 Running full notebook build + sync test"
	@make setup-test-block >/dev/null
	@make build
	@make edit-notebook
	@make reverse-sync
	@make verify-sync

build:
	@echo "📦 Building notebook from $(CONFIG)"
	../$(SCRIPT_DIR)/build.sh notebooks/differential_analysis/notebook_generator/tests/$(CONFIG)

setup-test-block:
	@mkdir -p blocks/test/
	@mkdir -p configs/
	@mkdir -p notebooks/
	@echo "{" > $(BLOCK)
	@echo "  \"cells\": [" >> $(BLOCK)
	@echo "   {" >> $(BLOCK)
	@echo "    \"cell_type\": \"markdown\"," >> $(BLOCK)
	@echo "    \"metadata\": { }," >> $(BLOCK)
	@echo "    \"source\": [" >> $(BLOCK)
	@echo "     \"Goal\"" >> $(BLOCK)
	@echo "    ]\n   }\n  ]\n }\n" >> $(BLOCK)
	@echo "block_root: notebooks/differential_analysis/notebook_generator/tests/blocks/" > $(CONFIG)
	@echo "projects:" >> $(CONFIG)
	@echo "  - name: TestNotebook" >> $(CONFIG)
	@echo "    blocks:" >> $(CONFIG)
	@echo "    - test/01_hello.json" >> $(CONFIG)
	git add $(BLOCK) && git commit -m "🧪 Add test block for sync" >/dev/null || true

verify-sync:
	@echo "🔍 Verifying that block was updated"
	@grep -q "Way" $(BLOCK) \
	&& echo "✅ Test passed: Block updated from notebook" \
	|| (echo "❌ Test failed: Block was not updated from notebook"; exit 1)

clean:
	@rm -f $(BLOCK) $(NOTEBOOK) $(CONFIG)
	@rm -rf ../blocks/test ../notebooks/TestNotebook.ipynb ../configs/test.yaml