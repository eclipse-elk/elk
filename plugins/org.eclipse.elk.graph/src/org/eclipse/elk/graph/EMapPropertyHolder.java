/**
 * Copyright (c) 2016 Kiel University and others.
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
package org.eclipse.elk.graph;

import org.eclipse.elk.graph.properties.IProperty;
import org.eclipse.elk.graph.properties.IPropertyHolder;
import org.eclipse.emf.common.util.EMap;

import org.eclipse.emf.ecore.EObject;

/**
 * <!-- begin-user-doc -->
 * A representation of the model object '<em><b>EMap Property Holder</b></em>'.
 * <!-- end-user-doc -->
 *
 * <!-- begin-model-doc -->
 * A property holder implementation based on {@link org.eclipse.emf.common.util.EMap} which can be used in Ecore models.
 * 
 * <p>This property holder implementation currently has two ways for saving properties: a map of properties as well as a map of <em>persistent entries</em>. Persistent entries are String-String pairs containing String representations of properties. When a graph is serialized, it is the persistent entries that get serialized, not the properties themselves. This has two implications. First, to save a graph, one has to call {@link #makePersistent()} first. Second, after loading a graph, one of the methods in {@link org.eclipse.elk.core.util.GraphDataUtil} needs to be called to turn persistent entries back into usable properties.</p>
 * <!-- end-model-doc -->
 *
 * <p>
 * The following features are supported:
 * </p>
 * <ul>
 *   <li>{@link org.eclipse.elk.graph.EMapPropertyHolder#getProperties <em>Properties</em>}</li>
 * </ul>
 *
 * @see org.eclipse.elk.graph.ElkGraphPackage#getEMapPropertyHolder()
 * @model abstract="true" superTypes="org.eclipse.elk.graph.IPropertyHolder"
 * @generated
 */
public interface EMapPropertyHolder extends EObject, IPropertyHolder {
    /**
     * Returns the value of the '<em><b>Properties</b></em>' map.
     * The key is of type {@link org.eclipse.elk.graph.properties.IProperty<?>},
     * and the value is of type {@link java.lang.Object},
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * <!-- begin-model-doc -->
     * Map of properties configured for this property holder.
     * <!-- end-model-doc -->
     * @return the value of the '<em>Properties</em>' map.
     * @see org.eclipse.elk.graph.ElkGraphPackage#getEMapPropertyHolder_Properties()
     * @model mapType="org.eclipse.elk.graph.ElkPropertyToValueMapEntry&lt;org.eclipse.elk.graph.IProperty&lt;?&gt;, org.eclipse.elk.graph.PropertyValue&gt;"
     * @generated
     */
    EMap<IProperty<?>, Object> getProperties();

} // EMapPropertyHolder
