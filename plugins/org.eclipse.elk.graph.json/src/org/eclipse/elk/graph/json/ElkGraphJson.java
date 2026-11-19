/**
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.

 * This Source Code may also be made available under the following Secondary
 * Licenses when the conditions for such availability set forth in the Eclipse
 * Public License v. 2.0 are satisfied: GPL-3.0 which is available at
 * https://www.gnu.org/licenses/gpl-3.0-standalone.html.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 */
package org.eclipse.elk.graph.json;

import java.io.IOException;
import java.io.StringReader;

import org.eclipse.elk.core.util.Maybe;
import org.eclipse.elk.graph.ElkNode;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonElement;
import com.google.gson.JsonIOException;
import com.google.gson.JsonObject;
import com.google.gson.TypeAdapter;
import com.google.gson.stream.JsonReader;

/**
 * Utility methods to import and export the ELK Graph JSON Format.
 */
public final class ElkGraphJson {

    private ElkGraphJson() {
    }

    /**
     * Initializes an importer for the passed json string. The import is finished using the
     * {@link ImportBuilder#toElk()} method.
     * 
     * @param graph
     *            the graph to import.
     * @return a builder instance that can be further configured.
     */
    public static ImportBuilder forGraph(final String graph) {
        ImportBuilder ib = new ImportBuilder();
        ib.graph = graph;
        return ib;
    }

    /**
     * Initializes an importer for the passed json object. The import is finished using the
     * {@link ImportBuilder#toElk()} method.
     * 
     * @param graph
     *            the graph to import.
     * @return a builder instance that can be further configured.
     * 
     * @see ImportBuilder#rememberImporter(Maybe)
     */
    public static ImportBuilder forGraph(final JsonObject graph) {
        ImportBuilder ib = new ImportBuilder();
        ib.jsonGraph = graph;
        return ib;
    }
    
    /**
     * Builder for importing.
     */
    public static final class ImportBuilder {
        
        private JsonObject jsonGraph;
        private String graph;
        private Maybe<JsonImporter> importerMaybe;
        /** See {@link JsonReader#setLenient(boolean)} for details. */
        private boolean lenient = true;
        
        
        /**
         * In case a later application of layout information is desired, pass a {@link Maybe} instance to this method.
         * The instance is populated with the used {@link JsonImporter}, which can be used to
         * {@link JsonImporter#transferLayout(ElkNode) transfer} the layout later on.
         * 
         * @param maybe
         *            an empty {@link Maybe}.
         * @return the builder.
         */
        public ImportBuilder rememberImporter(final Maybe<JsonImporter> maybe) {
            this.importerMaybe = maybe;
            return this;
        }

        /**
         * If <code>true</code>, be more relaxed when parsing the JSON, e.g. allow property keys without quotes. By
         * default <code>true</code>. For details see {@link JsonReader#setLenient(boolean)}.
         */
        public ImportBuilder lenient(final boolean beLenient) {
            this.lenient = beLenient;
            return this;
        }

        /**
         * Part of a GSON workaround as GSON is by default always 'lenient'. See
         * https://github.com/google/gson/issues/1208.
         */
        private static final TypeAdapter<JsonElement> GSON_ELEMENT_ADAPTER = new Gson().getAdapter(JsonElement.class);

        /**
         * Perform the actual import and return the resulting ELK Graph.
         * 
         * @return the root node of the imported ELK Graph.
         */
        public ElkNode toElk() {
            if (jsonGraph == null) {
                // Due to a GSON workaround the following lines are a bit more complicated that they have to be.
                // See the javadoc comment of GSON_ELEMENT_ADAPTER for details.
                JsonReader reader = new JsonReader(new StringReader(graph));
                reader.setLenient(this.lenient);
                try {
                    JsonElement json = GSON_ELEMENT_ADAPTER.read(reader);
                    if (!(json instanceof JsonObject)) {
                        throw new JsonImportException("Top-level element of the graph must be a json object.");
                    }
                    jsonGraph = json.getAsJsonObject();
                } catch (IOException e) {
                    throw new JsonIOException(e);
                }
            }

            JsonImporter importer = new JsonImporter();
            ElkNode elkGraph = importer.transform(jsonGraph);

            if (importerMaybe != null) {
                importerMaybe.set(importer);
            }

            return elkGraph;
        }
    }

    /**
     * Initializes the export of the passed graph. The export is finished using the {@link ExportBuilder#toJson()}
     * method.
     * 
     * @param graph
     *            the graph to export.
     * @return a builder instance that can be further configured.
     */
    public static ExportBuilder forGraph(final ElkNode graph) {
        ExportBuilder eb = new ExportBuilder();
        eb.graph = graph;
        return eb;
    }

    /**
     * Builder for exporting.
     */
    public static final class ExportBuilder {
        private ElkNode graph;
        private boolean prettyPrint = false;
        private boolean omitZeroPosition = true;
        private boolean omitZeroDimension = true;
        private boolean omitLayoutInformation = false;
        private boolean shortLayoutOptionKeys = true;
        private boolean omitUnknownLayoutOptions = true;

        /**
         * Pretty print the resulting json string, e.g. by inserting line breaks.
         */
        public ExportBuilder prettyPrint(final boolean pretty) {
            this.prettyPrint = pretty;
            return this;
        }

        /**
         * Omit any superfluous parts of the layout option keys, e.g. 'org.eclipse'.
         */
        public ExportBuilder shortLayoutOptionKeys(final boolean shortKeys) {
            this.shortLayoutOptionKeys = shortKeys;
            return this;
        }

        /**
         * Omit the position of an element if both are zero.
         */
        public ExportBuilder omitZeroPositions(final boolean omitZeroPos) {
            this.omitZeroPosition = omitZeroPos;
            return this;
        }

        /**
         * Omit the width and height of elements if they are zero.
         */
        public ExportBuilder omitZeroDimension(final boolean omitZeroDim) {
            this.omitZeroDimension = omitZeroDim;
            return this;
        }

        /**
         * Omit any layout information, e.g. node positions.
         */
        public ExportBuilder omitLayout(final boolean omitLayout) {
            this.omitLayoutInformation = omitLayout;
            return this;
        }

        /**
         * Omit any unknown layout options.
         */
        public ExportBuilder omitUnknownLayoutOptions(final boolean omitUnknownOptions) {
            this.omitUnknownLayoutOptions = omitUnknownOptions;
            return this;
        }
        
        /**
         * Perform the export using the specified configuration.
         * 
         * @return the json string representation of the graph.
         */
        public String toJson() {
            // export the graph
            JsonExporter exporter = new JsonExporter();
            exporter.setOptions(omitZeroPosition, omitZeroDimension, omitLayoutInformation,
                    shortLayoutOptionKeys, omitUnknownLayoutOptions);
            JsonObject jsonGraph = exporter.export(graph);

            // configure the gson builder
            GsonBuilder builder = new GsonBuilder();
            // no need to perform html escaping in our case, hence deactivate it:
            builder.disableHtmlEscaping();
            if (prettyPrint) {
                builder.setPrettyPrinting();
            }
            Gson gson = builder.create();

            // return the json string
            String json = gson.toJson(jsonGraph);
            return json;
        }
    }

}
