/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.util.persistence;

import java.util.Map;

import org.eclipse.emf.common.util.URI;
import org.eclipse.emf.ecore.xmi.XMLHelper;
import org.eclipse.emf.ecore.xmi.XMLSave;
import org.eclipse.emf.ecore.xmi.impl.XMIResourceImpl;

/**
 * A specialized resource for ELK graphs registered through extension points by the plug-in. This implementation
 * customizes the way ELK graphs are serialized and deserialized in the following ways.
 * 
 * 
 * <h3>Serialization</h3>
 * 
 * <p>When saving a graph, normally all of its elements and properties would be saved in a String representation.
 * To be deserialized again, however, Object type properties must be associated with a data type class which
 * implements {@link org.eclipse.elk.core.util.IDataObject IDataObject}. Not all properties do that. If a
 * property does not, it doesn't make sense to serialize it since we won't be able to deserialize them again
 * anyway. This implementation causes such property values to not be serialized.</p>
 * 
 * <p>Instances of this class default to UTF-8 encoding.</p>
 * 
 * 
 * <h3>Deserialization</h3>
 * 
 * <p>Deserializing properties requires the
 * {@link org.eclipse.elk.core.data.LayoutMetaDataService LayoutMetaDataService} to obtain appropriate
 * {@link org.eclipse.elk.core.data.LayoutOptionData LayoutOptionData} objects to be able to parse option
 * values. This cannot be done in the graph plug-in with the graph factory's conversion methods, since the
 * graph plug-in does not have access to the meta data service. Thus, this implementation catches these cases
 * and handles them properly.</p>
 * 
 * <p>Note that the old graph structure featured a <em>persistent entries</em> map which was a String-String
 * representation of properties and properties values. Using this implementation, that persistent entries map
 * became obsolete.</p>
 * 
 * @see ElkGraphResourceFactory
 * @see ElkGraphXMISave
 * @see ElkGraphXMIHelper
 */
public class ElkGraphResource extends XMIResourceImpl {
    
    /**
     * Creates an instance of the resource.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @param uri the URI of the new resource.
     * @generated
     */
    public ElkGraphResource(URI uri) {
        super(uri);
        
        // Default to utf-8 encoding (#163)
        setEncoding("utf-8");
    }
    
    
    @Override
    protected XMLSave createXMLSave() {
        return new ElkGraphXMISave(createXMLHelper());
    }
    
    @Override
    protected XMLSave createXMLSave(Map<?, ?> options) {
        // We simply always return our custom XMISave subclass regardless of whether or not XMI is to be suppressed
        return createXMLSave();
    }
    
    @Override
    protected XMLHelper createXMLHelper() {
        return new ElkGraphXMIHelper(this);
    }
    
}
